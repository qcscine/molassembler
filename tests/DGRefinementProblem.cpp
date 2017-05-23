#define BOOST_TEST_MODULE RefinementProblemTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "CommonTrig.h"
#include "DistanceGeometry/generateConformation.h"
#include "DistanceGeometry/MetricMatrix.h"
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

/* TODO
 */

bool isApprox(
  const dlib::matrix<double, 0, 1>& a,
  const dlib::matrix<double, 0, 1>& b,
  const double& epsilon
) {
  return TemplateMagic::all_of(
    TemplateMagic::zipMapAlternate(
      a,
      b,
      [&](const auto& i, const auto& j) -> bool {
        return std::fabs(
          std::fabs(i) - std::fabs(j)
        ) < epsilon;
      }
    )
  );
}

BOOST_AUTO_TEST_CASE( cppoptlibGradientCorrectnessCheck ) {
  using Vector = dlib::matrix<double, 0, 1>;

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

    // Transfer to dlib
    const Vector dlibPositions = dlib::mat(
      Eigen::VectorXd (
        Eigen::Map<Eigen::VectorXd>(
          embedded.data(),
          embedded.cols() * embedded.rows()
        )
      )
    );

    auto chiralityConstraints = TemplateMagic::map(
      DGInfo.chiralityConstraintPrototypes,
      [&DGInfo](
        const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
      ) -> ChiralityConstraint {
        return detail::propagate(
          DGInfo.distanceBounds,
          prototype
        );
      }
    );

    errfValue<false> valueFunctor {
      DGInfo.distanceBounds,
      chiralityConstraints
    };

    errfGradient<false> gradientFunctor {
      DGInfo.distanceBounds,
      chiralityConstraints
    };

    // Finite difference is calculated to 1e-7 precision
    Vector finiteDifferenceGradient = dlib::derivative(valueFunctor)(dlibPositions);
    Vector gradient = gradientFunctor(dlibPositions);

    bool passes = TemplateMagic::all_of(
      TemplateMagic::zipMapAlternate(
        gradient,
        finiteDifferenceGradient,
        [](const auto& a, const auto& b) -> bool {
          return std::fabs(
            std::fabs(a) - std::fabs(b)
          ) < 1e-5;
        }
      )
    );

    BOOST_CHECK_MESSAGE(
      passes,
      "Gradient does not equal numerical differentiation of value for "
        << Symmetry::name(symmetryName)
        << " error function without compression! Difference norm: "
        << dlib::length(gradient - finiteDifferenceGradient)
        << "\n positions: " << dlib::trans(dlibPositions)
        << "\n  gradient: " << dlib::trans(gradient)
        << "\n  finite diff grad: " << dlib::trans(finiteDifferenceGradient)
        << "\n  diff: " << dlib::trans(gradient - finiteDifferenceGradient)
    );

    errfValue<true> compressingValueFunctor {
      DGInfo.distanceBounds,
      chiralityConstraints
    };

    errfGradient<true> compressingGradientFunctor {
      DGInfo.distanceBounds,
      chiralityConstraints
    };

    Vector compressedFiniteDifferenceGradient = dlib::derivative(compressingValueFunctor)(dlibPositions);
    Vector compressedGradient = compressingGradientFunctor(dlibPositions);

    bool compressedPasses = TemplateMagic::all_of(
      TemplateMagic::zipMapAlternate(
        compressedGradient,
        compressedFiniteDifferenceGradient,
        [](const auto& a, const auto& b) -> bool {
          return std::fabs(
            std::fabs(a) - std::fabs(b)
          ) < 1e-5;
        }
      )
    );

    BOOST_CHECK_MESSAGE(
      compressedPasses,
      "Gradient does not equal numerical differentiation of value for "
        << Symmetry::name(symmetryName)
        << " compressed error function! Difference norm: "
        << dlib::length(compressedGradient - compressedFiniteDifferenceGradient)
        << "\n positions: " << dlib::trans(dlibPositions)
        << "\n  compressed gradient: " << dlib::trans(compressedGradient)
        << "\n  compressed finite diff grad: " << dlib::trans(compressedFiniteDifferenceGradient)
        << "\n  diff: " << dlib::trans(compressedGradient - compressedFiniteDifferenceGradient)
    );
  }
}

BOOST_AUTO_TEST_CASE( valueComponentsAreRotTransInvariant ) {
  for(const auto& symmetryName : Symmetry::allNames) {
    auto molecule = DGDBM::asymmetricMolecule(symmetryName);

    const auto DGData = gatherDGInformation(molecule);
    const auto distancesMatrix = DGData.distanceBounds.generateDistanceMatrix();
    auto chiralityConstraints = TemplateMagic::map(
      DGData.chiralityConstraintPrototypes,
      [&DGData](
        const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
      ) -> ChiralityConstraint {
        return detail::propagate(
          DGData.distanceBounds,
          prototype
        );
      }
    );

    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Transfer to dlib and vectorize
    errfValue<true>::Vector referencePositions = dlib::mat(
      static_cast<Eigen::VectorXd>(
          Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      )
    );

    errfValue<true> valueFunctor {
      DGData.distanceBounds,
      chiralityConstraints
    };

    // get a value
    const double referenceDistanceError = valueFunctor.distanceError(referencePositions);
    const double referenceChiralError = valueFunctor.chiralError(referencePositions);
    const double referenceExtraDimError = valueFunctor.extraDimensionError(referencePositions);

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

      // Transfer to dlib and vectorize
      errfValue<true>::Vector vectorizedPositions = dlib::mat(
        static_cast<Eigen::VectorXd>(
          Eigen::Map<Eigen::VectorXd>(
            embeddedCopy.data(),
            embeddedCopy.cols() * embeddedCopy.rows()
          )
        )
      );

      const double transformedDistanceError = valueFunctor.distanceError(vectorizedPositions);
      const double transformedChiralError = valueFunctor.chiralError(vectorizedPositions);
      const double transformedExtraDimError = valueFunctor.extraDimensionError(vectorizedPositions);

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
  using DlibVector = dlib::matrix<double, 0, 1>;

  for(const auto& symmetryName : Symmetry::allNames) {
    auto molecule = DGDBM::asymmetricMolecule(symmetryName);

    const auto DGData = gatherDGInformation(molecule);
    const auto distancesMatrix = DGData.distanceBounds.generateDistanceMatrix();
    auto chiralityConstraints = TemplateMagic::map(
      DGData.chiralityConstraintPrototypes,
      [&DGData](
        const Stereocenters::Stereocenter::ChiralityConstraintPrototype& prototype
      ) -> ChiralityConstraint {
        return detail::propagate(
          DGData.distanceBounds,
          prototype
        );
      }
    );

    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Vectorize the positions for use with cppoptlib
    errfValue<true>::Vector referencePositions = dlib::mat(
      static_cast<Eigen::VectorXd>(
          Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      )
    );

    errfGradient<true> gradientFunctor {
      DGData.distanceBounds,
      chiralityConstraints
    };

    const unsigned N = referencePositions.size() / 4;
    assert(N > 0);

    std::vector<DlibVector> referenceGradients;

    referenceGradients.emplace_back(gradientFunctor.referenceA(referencePositions));
    referenceGradients.emplace_back(gradientFunctor.referenceB(referencePositions));
    referenceGradients.emplace_back(gradientFunctor.referenceC(referencePositions));
    referenceGradients.emplace_back(gradientFunctor.referenceD(referencePositions));
    referenceGradients.emplace_back(gradientFunctor.referenceE(referencePositions));

    for(const auto& referenceGradient: referenceGradients) {
      assert(referenceGradient.size() == 4 * N);
    }

    /* Check that sum of reference implementations is equal to optimized
     * implementation 
     */
    DlibVector optimizedGradient = gradientFunctor(referencePositions);
    DlibVector referenceGradientsSum = (
      referenceGradients.at(0)
      + referenceGradients.at(1)
      + referenceGradients.at(2)
      + referenceGradients.at(3)
      + referenceGradients.at(4)
    );

    BOOST_CHECK_MESSAGE(
      isApprox(
        optimizedGradient,
        referenceGradientsSum,
        1e-10
      ),
      "Optimized gradient implementation is not equal to sum of reference gradients to 1e-10! Difference norm: "
      << dlib::length(referenceGradientsSum - optimizedGradient)
    );

    for(unsigned testNum = 0; testNum < 100; testNum++) {
      // Transformations
      const Eigen::Matrix3d rotationMatrix = Eigen::Quaterniond::UnitRandom().toRotationMatrix();

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
      errfValue<true>::Vector rotatedPositions = dlib::mat(
        static_cast<Eigen::VectorXd>(
          Eigen::Map<Eigen::VectorXd>(
            rotatedRectangularPositions.data(),
            rotatedRectangularPositions.cols() * rotatedRectangularPositions.rows()
          )
        )
      );

      errfValue<true>::Vector translatedPositions = dlib::mat(
        static_cast<Eigen::VectorXd>(
          Eigen::Map<Eigen::VectorXd>(
            translatedRectangularPositions.data(),
            translatedRectangularPositions.cols() * rotatedRectangularPositions.rows()
          )
        )
      );

      std::vector<DlibVector> rotatedGradients, translatedGradients;

      rotatedGradients.emplace_back(gradientFunctor.referenceA(rotatedPositions));
      rotatedGradients.emplace_back(gradientFunctor.referenceB(rotatedPositions));
      rotatedGradients.emplace_back(gradientFunctor.referenceC(rotatedPositions));
      rotatedGradients.emplace_back(gradientFunctor.referenceD(rotatedPositions));
      rotatedGradients.emplace_back(gradientFunctor.referenceE(rotatedPositions));

      translatedGradients.emplace_back(gradientFunctor.referenceA(translatedPositions));
      translatedGradients.emplace_back(gradientFunctor.referenceB(translatedPositions));
      translatedGradients.emplace_back(gradientFunctor.referenceC(translatedPositions));
      translatedGradients.emplace_back(gradientFunctor.referenceD(translatedPositions));
      translatedGradients.emplace_back(gradientFunctor.referenceE(translatedPositions));

      // Transform the rotated gradients
      auto rotatedReferenceGradients = referenceGradients;

      // Transfer rotation matrix to dlib
      dlib::matrix<double, 3, 3> dlibRotationMatrix = dlib::mat(rotationMatrix);

      for(auto& rotatedReferenceGradient : rotatedReferenceGradients) {
        for(unsigned i = 0; i < N; i++) {
          dlib::set_rowm(
            rotatedReferenceGradient,
            dlib::range(4 * i, 4 * i + 2)
          ) = dlibRotationMatrix * dlib::rowm(  
            rotatedReferenceGradient,
            dlib::range(4 * i, 4 * i + 2)
          );
        }
      }

      // Compare
      for(unsigned i = 0; i < 5; i++) {
        BOOST_CHECK_MESSAGE(
          isApprox(
            rotatedReferenceGradients[i],
            rotatedGradients[i],
            1e-10
          ),
          "Gradient component " << static_cast<char>('A' + i) 
            << " is not rotationally invariant! Difference norm: "
            << dlib::length(rotatedReferenceGradients[i] - rotatedGradients[i])
        );

        BOOST_CHECK_MESSAGE(
          isApprox(
            referenceGradients[i],
            translatedGradients[i],
            1e-10
          ),
          "Gradient component " << static_cast<char>('A' + i) 
            << " is not translationally invariant! Difference norm: "
            << dlib::length(referenceGradients[i] - translatedGradients[i])
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

  const double maximumErrorThreshold = 1e-5;

  // Open output file for R plots
  std::string filename = "DGRefinementProblem-symmetric-ensemble-errors.csv";
  std::ofstream outFile(filename);
  outFile << std::setprecision(4) << std::fixed;

  for(const auto& enumPair : enumerate(Symmetry::allNames)) {
    const auto& symmetryName = enumPair.value;
    const auto& symmetryIndex = enumPair.index;

    std::cout << Symmetry::name(symmetryName) << std::endl;

    auto molecule = DGDBM::asymmetricMolecule(symmetryName);
    
    auto DGResult = DistanceGeometry::detail::debugDistanceGeometry(
      molecule,
      100,
      MetrizationOption::full,
      false,
      BFSConstraintCollector::DistanceMethod::Uniform
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

    // The average error of the ensemble should be below 1e-5 (already achieved)
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
