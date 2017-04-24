#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/MetricMatrix.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "cppoptlib/meta.h"

#include "BoundsFromSymmetry.h"
#include "IO.h"

#include "AnalysisHelpers.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

/* NOTES to current state of optimized output structures
 *
 * - Initial problem is gone, was an issue with the generation of the test
 *   matrices
 * - Found one bug in the refinement problem (a missing square), heavily
 *   decreases optimization time
 * - Found another bug that caused generated structures to be considerably
 *   expanded, although well reflective of the overall symmetries (a missing
 *   sqrt in the embedding procedure).
 * - Now high-symmetry geometries are maybe somewhat over-constrained, leading 
 *   to messy structures. Avenues to try:
 *
 *   - Add more variance to 1-3 distance bounds
 *     -> did nothing
 *
 *   - Add more sanity tests to the various symmetries' angle functions
 *     -> Found mistakes in the angle functions of pent-bipy and sq-antiprism, 
 *        fixing the issue
 *        
 * - Refinement cannot reach a global minimum in almost all cases. Probably due
 *   to triangle inequality bound violations. Fix by introducing metrization.
 *
 * - Chirality constraints are now active. A variation of how they are defined 
 *   exists in the symmetry_information library (if nothing else helps to fix 
 *   the appearing bugs). There was a mistake in the Cayley-Menger determinant 
 *   calculation to determine the target volumes, that's fixed. There's now a
 *   weird bug where gradients with chirality constraints are quite unstable.
 */

int main() {

  auto embeddingOption = EmbeddingOption::fourDimensional;

  // Make a problem solver
  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGSolver;

  // Set stop criteria
  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  //stopCriteria.iterations = 1000;
  stopCriteria.fDelta = 1e-5;
  DGSolver.setStopCriteria(stopCriteria);

  const unsigned nStructures = 100;

  for(const auto& symmetryName : Symmetry::allNames) {
    // Make a space-free string from the name
    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    // Generate distance bounds
    auto simpleMol = DGDBM::symmetricMolecule(symmetryName);

    auto DGInfo = DistanceGeometry::gatherDGInformation(simpleMol);
    const auto& distanceBoundsMatrix = DGInfo.distanceBounds;

    /*std::cout << "Sample distances matrix for symmetry '" 
      << Symmetry::name(symmetryName) << std::endl
      << distanceBoundsMatrix.generateDistanceMatrix(
        MetrizationOption::off
      ) << std::endl;*/

    std::vector<double> refinedErrorValues;

    for(unsigned structNum = 0; structNum < nStructures; structNum++) {

      auto distancesMatrix = distanceBoundsMatrix.generateDistanceMatrix(); // with metrization!

      // Calculate metric matrix from selected distances
      MetricMatrix metricMatrix(distancesMatrix);

      // Embed
      auto embeddedPositions = metricMatrix.embed(embeddingOption);

      // Vectorize
      Eigen::VectorXd vectorizedPositions(
        Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      );

      // Write the unoptimized result to a file
      writeMOLFile(
        simpleMol,
        false,
        spaceFreeName,
        structNum,
        vectorizedPositions,
        embeddingOption
      );

      // Get constraints
      auto chiralityConstraints = TemplateMagic::map(
        DGInfo.chiralityConstraintPrototypes,
        /* Partial application with distances matrix so we have a unary function
         * to perform the mapping with
         */
        detail::makePropagator(
          [&distancesMatrix](const AtomIndexType& i, const AtomIndexType& j) {
            return distancesMatrix(
              std::min(i, j),
              std::max(i, j)
            );
          }
        )
      );

      // Create the RefinementProblem
      DGRefinementProblem<double> problem(
        chiralityConstraints,
        distanceBoundsMatrix
      );

      // Run the minimization
      DGSolver.minimize(problem, vectorizedPositions);

      // Save end value in Problem
      refinedErrorValues.emplace_back(
        problem.value(
          vectorizedPositions
        )
      );

      // Write the result to a file
      writeMOLFile(
        simpleMol,
        true,
        spaceFreeName,
        structNum,
        vectorizedPositions,
        embeddingOption
      );
    }

    // Write refined error values to file
    writeErrorValues(
      spaceFreeName,
      refinedErrorValues
    );
  }
}
