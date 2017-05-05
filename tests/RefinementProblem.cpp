#define BOOST_TEST_MODULE RefinementProblemTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/MetricMatrix.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "cppoptlib/meta.h"
#include "DistanceGeometry/RefinementProblem.h"
#include "template_magic/Enumerate.h"
#include "template_magic/Random.h"
#include "AnalysisHelpers.h"

#include <Eigen/Geometry>

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
  
  // Generate a wide array of RefinementProblems and check their gradients

  for(const auto& symmetryName: Symmetry::allNames) {
    auto molecule = DGDBM::asymmetricMolecule(symmetryName);
    auto DGInfo = gatherDGInformation(molecule);

    auto distances = DGInfo.distanceBounds.generateDistanceMatrix(
      MetrizationOption::full
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

    RefinementProblem problem(
      chiralityConstraints,
      DGInfo.distanceBounds
    );

    BOOST_CHECK(
      problem.checkGradient(vectorizedPositions)
    );

    // TODO minimize, then check again?
  }
}

BOOST_AUTO_TEST_CASE( valueComponentsAreRotTransInvariant ) {
  for(const auto& symmetryName : Symmetry::allNames) {
    auto molecule = DGDBM::symmetricMolecule(symmetryName);

    const auto DGData = gatherDGInformation(molecule);
    const auto distancesMatrix = DGData.distanceBounds.generateDistanceMatrix();
    auto chiralityConstraints = TemplateMagic::map(
      DGData.chiralityConstraintPrototypes,
      detail::makePropagator(
        [&distancesMatrix](
          const AtomIndexType& i,
          const AtomIndexType& j
        ) -> double {
          return distancesMatrix(
            std::min(i, j),
            std::max(i, j)
          );
        }
      )
    );

    RefinementProblem problem(
      chiralityConstraints,
      DGData.distanceBounds
    );
    
    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Vectorize the positions for use with cppoptlib
    Eigen::VectorXd referencePositions {
      Eigen::Map<Eigen::VectorXd>(
        embeddedPositions.data(),
        embeddedPositions.cols() * embeddedPositions.rows()
      )
    };

    // get a value
    const double referenceDistanceError = problem.distanceError(referencePositions);
    const double referenceChiralError = problem.chiralError(referencePositions);
    const double referenceExtraDimError = problem.extraDimensionError(referencePositions);

    bool distanceErrorRotTransInvariant = true;
    bool chiralErrorRotTransInvariant = true;
    bool extraDimErrorRotTransInvariant = true;

    for(unsigned testNum = 0; testNum < 100; testNum++) {
      const Eigen::Matrix3d rotationMatrix = Eigen::Quaterniond::UnitRandom().toRotationMatrix();

      const Eigen::Vector3d translation = Eigen::Vector3d::Random();
      auto embeddedCopy = embeddedPositions;

      // Transform the embedded coordinates by rotation and translation
      for(unsigned i = 0; i < embeddedCopy.cols(); i++) {
        embeddedCopy.template block<3, 1>(0, i) = (
          rotationMatrix * embeddedCopy.template block<3, 1>(0, i)
        );
        embeddedCopy.template block<3, 1>(0, i) += translation;
      }

      // Vectorize the positions
      Eigen::VectorXd vectorizedPositions {
        Eigen::Map<Eigen::VectorXd>(
          embeddedCopy.data(),
          embeddedCopy.cols() * embeddedCopy.rows()
        )
      };

      const double transformedDistanceError = problem.distanceError(vectorizedPositions);
      const double transformedChiralError = problem.chiralError(vectorizedPositions);
      const double transformedExtraDimError = problem.extraDimensionError(vectorizedPositions);

      if(
        distanceErrorRotTransInvariant
        && std::fabs(transformedDistanceError - referenceDistanceError) > 1e-10
      ) {
        std::cout << "RefinementProblem distance error is not 3D rot-trans "
          << "invariant for " << Symmetry::name(symmetryName) << "." << std::endl;
        
        distanceErrorRotTransInvariant = false;
      }

      if(
        chiralErrorRotTransInvariant
        && std::fabs(transformedChiralError - referenceChiralError) > 1e-10
      ) {
        std::cout << "RefinementProblem chiral error is not 3D rot-trans "
          << "invariant for " << Symmetry::name(symmetryName) << "." << std::endl;
        chiralErrorRotTransInvariant = false;
      }

      if(
        extraDimErrorRotTransInvariant
        && std::fabs(transformedExtraDimError - referenceExtraDimError) > 1e-10
      ) {
        std::cout << "RefinementProblem extra dim error is not 3D rot-trans "
          << "invariant for " << Symmetry::name(symmetryName)
          << ". The fourth dimension is untouched in 3D rotation/translation!"
          << std::endl;
        extraDimErrorRotTransInvariant = false;
      }
    }

    BOOST_CHECK(distanceErrorRotTransInvariant);
    BOOST_CHECK(chiralErrorRotTransInvariant);
    BOOST_CHECK(extraDimErrorRotTransInvariant);
  }
}

BOOST_AUTO_TEST_CASE( gradientComponentsAreRotAndTransInvariant) {
  for(const auto& symmetryName : Symmetry::allNames) {
    auto molecule = DGDBM::symmetricMolecule(symmetryName);

    const auto DGData = gatherDGInformation(molecule);
    const auto distancesMatrix = DGData.distanceBounds.generateDistanceMatrix();
    auto chiralityConstraints = TemplateMagic::map(
      DGData.chiralityConstraintPrototypes,
      detail::makePropagator(
        [&distancesMatrix](
          const AtomIndexType& i,
          const AtomIndexType& j
        ) -> double {
          return distancesMatrix(
            std::min(i, j),
            std::max(i, j)
          );
        }
      )
    );

    RefinementProblem problem(
      chiralityConstraints,
      DGData.distanceBounds
    );

    // Set compress to avoid any state change issues
    problem.compress = true;

    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Vectorize the positions for use with cppoptlib
    const Eigen::VectorXd referencePositions {
      Eigen::Map<Eigen::VectorXd>(
        embeddedPositions.data(),
        embeddedPositions.cols() * embeddedPositions.rows()
      )
    };

    const Eigen::VectorXd emptyGradient = Eigen::VectorXd::Zero(referencePositions.size());

    const unsigned N = referencePositions.size() / 4;

    std::vector<Eigen::VectorXd> referenceGradients {4, emptyGradient};

    problem.gradientA(referencePositions, referenceGradients[0]);
    problem.gradientB(referencePositions, referenceGradients[1]);
    problem.gradientC(referencePositions, referenceGradients[2]);
    problem.gradientD(referencePositions, referenceGradients[3]);

    // TODO remove: Is split implemented the same as gradient? 
    auto fullGradient = emptyGradient;
    problem.gradient(referencePositions, fullGradient);
    BOOST_CHECK_MESSAGE(
      fullGradient.isApprox(
        referenceGradients[0]
        + referenceGradients[1]
        + referenceGradients[2]
        + referenceGradients[3],
        1e-10
      ),
      "Full gradient implementation is not equal to sum of split gradients to 1e-10!"
    );

    for(unsigned testNum = 0; testNum < 1; testNum++) {
      // Transformations
      const Eigen::Matrix3d rotationMatrix = 
         Eigen::Quaterniond::UnitRandom().toRotationMatrix();

      const Eigen::Vector3d translation = Eigen::Vector3d::Random();
      Eigen::MatrixXd rotatedRectangularPositions = embeddedPositions;
      Eigen::MatrixXd translatedRectangularPositions = embeddedPositions;

      // Transform the embedded coordinates by rotation and translation
      for(unsigned i = 0; i < rotatedRectangularPositions.cols(); i++) {
        rotatedRectangularPositions.template block<3, 1>(0, i) = (
          rotationMatrix * rotatedRectangularPositions.template block<3, 1>(0, i)
        );
        translatedRectangularPositions.template block<3, 1>(0, i) += translation;
      }

      // Vectorize the positions
      const Eigen::VectorXd rotatedPositions {
        Eigen::Map<Eigen::VectorXd>(
          rotatedRectangularPositions.data(),
          rotatedRectangularPositions.cols() * rotatedRectangularPositions.rows()
        )
      };

      const Eigen::VectorXd translatedPositions {
        Eigen::Map<Eigen::VectorXd>(
          translatedRectangularPositions.data(),
          translatedRectangularPositions.cols() * rotatedRectangularPositions.rows()
        )
      };

      std::vector<Eigen::VectorXd> rotatedGradients {4, emptyGradient};

      problem.gradientA(rotatedPositions, rotatedGradients[0]);
      problem.gradientB(rotatedPositions, rotatedGradients[1]);
      problem.gradientC(rotatedPositions, rotatedGradients[2]);
      problem.gradientD(rotatedPositions, rotatedGradients[3]);

      std::vector<Eigen::VectorXd> translatedGradients {4, emptyGradient};

      problem.gradientA(translatedPositions, translatedGradients[0]);
      problem.gradientB(translatedPositions, translatedGradients[1]);
      problem.gradientC(translatedPositions, translatedGradients[2]);
      problem.gradientD(translatedPositions, translatedGradients[3]);

      // Transform the rotated gradients
      std::vector<Eigen::VectorXd> rotatedReferenceGradients = referenceGradients;

      for(auto& rotatedReferenceGradient : rotatedReferenceGradients) {
        for(unsigned i = 0; i < N; i++) {
          rotatedReferenceGradient.template segment<3>(4 * i) = (
            rotationMatrix * rotatedReferenceGradient.template segment<3>(4 * i)
          );
        }
      }

      // Compare
      for(unsigned i = 0; i < 4; i++) {
        BOOST_CHECK_MESSAGE(
          rotatedReferenceGradients[i].isApprox(rotatedGradients[i], 1e-10)
          || (rotatedReferenceGradients[i] - rotatedGradients[i]).norm() < 1e-10,
          "Gradient component " << static_cast<char>('A' + i) 
            << " is not rotationally invariant! Difference norm: "
            << (rotatedReferenceGradients[i] - rotatedGradients[i]).norm()
        );

        BOOST_CHECK_MESSAGE(
          referenceGradients[i].isApprox(translatedGradients[i], 1e-10)
          || (referenceGradients[i] - translatedGradients[i]).norm() < 1e-10,
          "Gradient component " << static_cast<char>('A' + i) 
            << " is not translationally invariant! Difference norm: "
            << (referenceGradients[i] - translatedGradients[i]).norm()
        );
      }
    }
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
      MetrizationOption::full
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
