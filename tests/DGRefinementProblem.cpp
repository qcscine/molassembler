#define BOOST_TEST_MODULE DGRefinementProblemTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/MetricMatrix.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "cppoptlib/meta.h"
#include "DistanceGeometry/DGRefinementProblem.h"
#include "template_magic/Enumerate.h"
#include "AnalysisHelpers.h"

#include "BoundsFromSymmetry.h"
#include "IO.h"

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

/* TODO
 * - Since we get a bug before we can even confirm whether some sample molecules
 *   give the correct gradients, it might be worth adding more asserts() as 
 *   checks in the distance matrix generating code in generateConformation as
 *   well to ensure that e.g. distance matrices adhere to triangle inequalities,
 *   and more if you can think of any
 *
 *   Progress: Bug is basically that the trig prismatic matrix is unsmoothed, 
 *   has some default 100 distances in there for some reason I cannot fathom. 
 *   Check it and add more checks.
 */

BOOST_AUTO_TEST_CASE( cppoptlibGradientCorrectnessCheck ) {
  Log::level = Log::Level::None;
  Log::particulars = {};
  
  // Generate a wide array of DGRefinementProblems and check their gradients

  for(const auto& symmetryName: Symmetry::allNames) {
    auto molecule = DGDBM::asymmetricMolecule(symmetryName);
    auto DGInfo = gatherDGInformation(molecule);

    auto distances = DGInfo.distanceBounds.generateDistanceMatrix(
      MetrizationOption::off
    );

    MetricMatrix metric(distances);

    auto embedded = metric.embed();

    Eigen::VectorXd vectorizedPositions(
      Eigen::Map<Eigen::VectorXd>(
        embedded.data(),
        embedded.cols() * embedded.rows()
      )
    );

    auto chiralityConstraints = TemplateMagic::map(
      DGInfo.chiralityConstraintPrototypes,
      detail::makePropagator(
        [&distances](const AtomIndexType& i, const AtomIndexType& j) {
          return distances(
            std::min(i, j),
            std::max(i, j)
          );
        }
      )
    );

    DGRefinementProblem<double> problem(
      chiralityConstraints,
      DGInfo.distanceBounds
    );

    BOOST_CHECK(
      problem.checkGradient(vectorizedPositions)
    );

    // TODO minimize, then check again?
  }
}

BOOST_AUTO_TEST_CASE( basicMoleculeDGWorksWell ) {
  Log::level = Log::Level::None;
  Log::particulars = {
    Log::Particulars::DGDebugInfo,
    Log::Particulars::DGRefinementChiralityNumericalDebugInfo
  };

  const double maximumErrorThreshold = 0.1;

  // Open output file for R plots
  std::string filename = "DGRefinementProblem-symmetric-ensemble-errors.csv";
  std::ofstream outFile(filename);
  outFile << std::setprecision(4) << std::fixed;

  for(const auto& enumPair : enumerate(Symmetry::allNames)) {
    const auto& symmetryName = enumPair.value;
    const auto& symmetryIndex = enumPair.index;

    std::cout << Symmetry::name(symmetryName) << std::endl;

    auto molecule = DGDBM::symmetricMolecule(symmetryName);
    
    auto DGResult = DistanceGeometry::detail::debugDistanceGeometry(
      molecule,
      100,
      MetrizationOption::off
    );

    // For something this simple, there really shouldn't be any failures
    BOOST_CHECK(DGResult.failures == 0);

    auto sumErrors = [](const detail::RefinementStepData& stepData) -> double {
      return (
        stepData.distanceError
        + stepData.chiralError
        + stepData.fourthDimError
      );
    };

    auto finalErrors = TemplateMagic::map(
      DGResult.refinements,
      [&](const detail::RefinementData& refinementData) -> double {
        return sumErrors(refinementData.steps.back());
      }
    );

    // The average error of the ensemble should be below 0.1 (already achieved)
    BOOST_CHECK(TemplateMagic::numeric::average(finalErrors) < maximumErrorThreshold);

    for(const auto& enumPair : enumerate(DGResult.refinements)) {
      const auto& refinementData = enumPair.value;

      const auto& finalError = sumErrors(refinementData.steps.back());

      // Write the final error to file along with the symmetry index
      outFile << symmetryIndex << "," << finalError << "\n";

      // Should it exceed the threshold, write the corresponding debug files
      if(finalError >= maximumErrorThreshold) {
        AnalysisHelpers::writeDGPOVandProgressFiles(
          molecule,
          symmetryName,
          enumPair.index,
          refinementData
        );
      }
    }

  }

  outFile.close();
}
