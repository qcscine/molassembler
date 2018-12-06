/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "Eigen/Geometry"
#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"
#include "temple/Adaptors/Zip.h"
#include "temple/Functional.h"
#include "temple/Random.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/constexpr/Numeric.h"

#include "molassembler/Detail/AnalysisHelpers.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"
#include "molassembler/DistanceGeometry/MetricMatrix.h"
#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/RefinementProblem.h"
#include "molassembler/IO.h"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std::string_literals;
using namespace molassembler;
using namespace molassembler::DistanceGeometry;

bool isApprox(
  const dlib::matrix<double, 0, 1>& a,
  const dlib::matrix<double, 0, 1>& b,
  const double epsilon
) {
  return temple::all_of(
    temple::adaptors::zip(a, b),
    [&](const auto i, const auto j) -> bool {
      return temple::floating::isCloseAbsolute(
        i,
        j,
        epsilon
      );
    }
  );
}

BOOST_AUTO_TEST_CASE( numericDifferentiationApproximatesGradients ) {
  using Vector = dlib::matrix<double, 0, 1>;

  Log::level = Log::Level::None;
  Log::particulars = {};

  // Generate several RefinementProblems and check their gradients
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("ez_stereocenters")
  ) {
    auto molecule = IO::read(currentFilePath.string());

    auto DGInfo = gatherDGInformation(molecule, DistanceGeometry::Configuration {});

    DistanceBoundsMatrix distanceBounds {
      molecule,
      DGInfo.bounds
    };

    auto distancesResult = distanceBounds.makeDistanceMatrix();
    if(!distancesResult) {
      BOOST_FAIL(distancesResult.error().message());
    }
    auto distances = distancesResult.value();

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

    dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      static_cast<Eigen::MatrixXd>(
        distanceBounds.access().cwiseProduct(distanceBounds.access())
      )
    );

    errfValue<false> valueFunctor {
      squaredBounds,
      DGInfo.chiralityConstraints,
      DGInfo.dihedralConstraints
    };

    errfGradient<false> gradientFunctor {
      squaredBounds,
      DGInfo.chiralityConstraints,
      DGInfo.dihedralConstraints
    };

    // Finite difference is calculated to 1e-7 precision
    Vector finiteDifferenceGradient = dlib::derivative(valueFunctor)(dlibPositions);
    Vector gradient = gradientFunctor(dlibPositions);

    bool passes = temple::all_of(
      temple::adaptors::zip(
        gradient,
        finiteDifferenceGradient
      ),
      [](const double a, const double b) -> bool {
        return temple::floating::detail::isCloseRelativeOrAbsolute(
          a,
          b,
          1e-5,
          1e-5
        );
      }
    );

    BOOST_CHECK_MESSAGE(
      passes,
      "Gradient does not equal numerical differentiation of value for "
        << currentFilePath.string()
        << " error function without compression! Difference norm: "
        << dlib::length(gradient - finiteDifferenceGradient)
        << "\n positions: " << dlib::trans(dlibPositions)
        << "\n  gradient: " << dlib::trans(gradient)
        << "\n  finite diff grad: " << dlib::trans(finiteDifferenceGradient)
        << "\n  diff: " << dlib::trans(gradient - finiteDifferenceGradient)
    );

    errfValue<true> compressingValueFunctor {
      squaredBounds,
      DGInfo.chiralityConstraints,
      DGInfo.dihedralConstraints
    };

    errfGradient<true> compressingGradientFunctor {
      squaredBounds,
      DGInfo.chiralityConstraints,
      DGInfo.dihedralConstraints
    };

    Vector compressedFiniteDifferenceGradient = dlib::derivative(compressingValueFunctor)(dlibPositions);
    Vector compressedGradient = compressingGradientFunctor(dlibPositions);

    bool compressedPasses = temple::all_of(
      temple::adaptors::zip(
        compressedGradient,
        compressedFiniteDifferenceGradient
      ),
      [](const double a, const double b) -> bool {
        return temple::floating::detail::isCloseRelativeOrAbsolute(
          a,
          b,
          1e-5,
          1e-5
        );
      }
    );

    BOOST_CHECK_MESSAGE(
      compressedPasses,
      "Gradient does not equal numerical differentiation of value for "
        << currentFilePath.string()
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
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("ez_stereocenters")
  ) {
    auto molecule = IO::read(currentFilePath.string());

    const auto DGData = gatherDGInformation(molecule, DistanceGeometry::Configuration {});

    DistanceBoundsMatrix distanceBounds {
      molecule,
      DGData.bounds
    };

    auto distancesMatrixResult = distanceBounds.makeDistanceMatrix();
    if(!distancesMatrixResult) {
      BOOST_FAIL(distancesMatrixResult.error().message());
    }
    auto distancesMatrix = distancesMatrixResult.value();

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

    dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      static_cast<Eigen::MatrixXd>(
        distanceBounds.access().cwiseProduct(distanceBounds.access())
      )
    );

    errfValue<true> valueFunctor {
      squaredBounds,
      DGData.chiralityConstraints,
      DGData.dihedralConstraints
    };

    // get a value
    const double referenceDistanceError = valueFunctor.distanceError(referencePositions);
    const double referenceChiralError = valueFunctor.chiralError(referencePositions);
    const double referenceDihedralError = valueFunctor.dihedralError(referencePositions);
    const double referenceExtraDimError = valueFunctor.extraDimensionError(referencePositions);

    bool distanceErrorRotTransInvariant = true;
    bool chiralErrorRotTransInvariant = true;
    bool dihedralErrorRotTransInvariant = true;
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
      const double transformedDihedralError = valueFunctor.dihedralError(vectorizedPositions);
      const double transformedExtraDimError = valueFunctor.extraDimensionError(vectorizedPositions);

      if(
        distanceErrorRotTransInvariant
        && std::fabs(transformedDistanceError - referenceDistanceError) > 1e-10
      ) {
        std::cout << "RefinementProblem distance error is not 3D rot-trans "
          << "invariant for " << currentFilePath.string() << "." << std::endl;

        distanceErrorRotTransInvariant = false;
      }

      if(
        chiralErrorRotTransInvariant
        && std::fabs(transformedChiralError - referenceChiralError) > 1e-10
      ) {
        std::cout << "RefinementProblem chiral error is not 3D rot-trans "
          << "invariant for " << currentFilePath.string() << "." << std::endl;
        chiralErrorRotTransInvariant = false;
      }

      if(
        dihedralErrorRotTransInvariant
        && std::fabs(transformedDihedralError - referenceDihedralError) > 1e-10
      ) {
        std::cout << "RefinementProblem dihedral error is not 3D rot-trans "
          << "invariant for " << currentFilePath.string() << "." << std::endl;
        dihedralErrorRotTransInvariant = false;
      }

      if(
        extraDimErrorRotTransInvariant
        && std::fabs(transformedExtraDimError - referenceExtraDimError) > 1e-10
      ) {
        std::cout << "RefinementProblem extra dim error is not 3D rot-trans "
          << "invariant for " << currentFilePath.string()
          << ". The fourth dimension is untouched in 3D rotation/translation!"
          << std::endl;
        extraDimErrorRotTransInvariant = false;
      }
    }

    BOOST_CHECK(distanceErrorRotTransInvariant);
    BOOST_CHECK(chiralErrorRotTransInvariant);
    BOOST_CHECK(dihedralErrorRotTransInvariant);
    BOOST_CHECK(extraDimErrorRotTransInvariant);
  }
}

BOOST_AUTO_TEST_CASE( gradientComponentsAreRotAndTransInvariant) {
  using DlibVector = dlib::matrix<double, 0, 1>;

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("ez_stereocenters")
  ) {
    auto molecule = IO::read(currentFilePath.string());

    const auto DGData = gatherDGInformation(molecule, DistanceGeometry::Configuration {});
    DistanceBoundsMatrix distanceBounds {
      molecule,
      DGData.bounds
    };

    auto distancesMatrixResult = distanceBounds.makeDistanceMatrix();
    if(!distancesMatrixResult) {
      BOOST_FAIL(distancesMatrixResult.error().message());
    }
    auto distancesMatrix = distancesMatrixResult.value();

    // Make a metric matrix from the distances matrix
    MetricMatrix metric(distancesMatrix);

    // Get a position matrix by embedding the metric matrix
    auto embeddedPositions = metric.embed();

    // Vectorize the positions for use with dlib
    errfValue<true>::Vector referencePositions = dlib::mat(
      static_cast<Eigen::VectorXd>(
          Eigen::Map<Eigen::VectorXd>(
          embeddedPositions.data(),
          embeddedPositions.cols() * embeddedPositions.rows()
        )
      )
    );

    dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      static_cast<Eigen::MatrixXd>(
        distanceBounds.access().cwiseProduct(distanceBounds.access())
      )
    );

    errfGradient<true> gradientFunctor {
      squaredBounds,
      DGData.chiralityConstraints,
      DGData.dihedralConstraints
    };

    const unsigned N = referencePositions.size() / 4;
    assert(N > 0);

    std::vector<DlibVector> referenceGradients;

    referenceGradients.emplace_back(gradientFunctor.referenceA(referencePositions));
    referenceGradients.emplace_back(gradientFunctor.referenceB(referencePositions));
    referenceGradients.emplace_back(gradientFunctor.referenceC(referencePositions));
    referenceGradients.emplace_back(gradientFunctor.referenceD(referencePositions));
    referenceGradients.emplace_back(gradientFunctor.referenceDihedral(referencePositions));

    assert(
      temple::all_of(
        referenceGradients,
        [&N](const auto& referenceGradient) -> bool {
          return referenceGradient.size() == 4 * N;
        }
      )
    );

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
      rotatedGradients.emplace_back(gradientFunctor.referenceDihedral(rotatedPositions));

      translatedGradients.emplace_back(gradientFunctor.referenceA(translatedPositions));
      translatedGradients.emplace_back(gradientFunctor.referenceB(translatedPositions));
      translatedGradients.emplace_back(gradientFunctor.referenceC(translatedPositions));
      translatedGradients.emplace_back(gradientFunctor.referenceD(translatedPositions));
      translatedGradients.emplace_back(gradientFunctor.referenceDihedral(translatedPositions));

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
      for(unsigned i = 0; i < 4; i++) {
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
