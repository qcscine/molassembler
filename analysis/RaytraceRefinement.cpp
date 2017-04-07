#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
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

int main() {

  // Make a problem solver
  cppoptlib::ConjugatedGradientDescentSolver<
    DGRefinementProblem<double>
  > DGSolver;

  // Set stop criteria
  cppoptlib::Criteria<double> stopCriteria = cppoptlib::Criteria<double>::defaults();
  stopCriteria.iterations = 1000;
  stopCriteria.fDelta = 1e-5;
  DGSolver.setStopCriteria(stopCriteria);

  const unsigned nStructures = 10;

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
    auto distanceBoundsMatrix = simpleMol.getDistanceBoundsMatrix();

    std::vector<double> refinedErrorValues;

    for(unsigned structNum = 0; structNum < nStructures; structNum++) {

      // Calculate metric matrix from selected distances
      MetricMatrix metricMatrix(
        distanceBoundsMatrix.generateDistanceMatrix() // with metrization!
      );

      // Embed
      auto embeddedPositions = metricMatrix.embed(EmbeddingOption::threeDimensional);

      // Vectorize
      Eigen::VectorXd vectorizedPositions(
        Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      );

      // Create the RefinementProblem
      DGRefinementProblem<double> problem(
        std::vector<ChiralityConstraint>({}), // no chirality constraints
        distanceBoundsMatrix
      );

      // Run the minimization
      auto stepResult = DGSolver.step(problem, vectorizedPositions);
      unsigned iterations = 1;
      writePOVRayFile(
        spaceFreeName,
        structNum,
        iterations,
        problem,
        simpleMol,
        vectorizedPositions,
        stepResult
      );

      while(stepResult.status == cppoptlib::Status::Continue) {
        stepResult = DGSolver.step(problem, vectorizedPositions);
        iterations += 1;
        writePOVRayFile(
          spaceFreeName,
          structNum,
          iterations,
          problem,
          simpleMol,
          vectorizedPositions,
          stepResult
        );
      }
    }
  }
}
