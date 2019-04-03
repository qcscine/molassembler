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
#include "temple/Functor.h"
#include "temple/Random.h"
#include "temple/Stringify.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/constexpr/Numeric.h"
#include "temple/constexpr/TupleType.h"

#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"
#include "molassembler/DistanceGeometry/MetricMatrix.h"
#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/DlibRefinement.h"
#include "molassembler/DistanceGeometry/EigenRefinement.h"
#include "molassembler/DistanceGeometry/EigenSIMDRefinement.h"
#include "molassembler/IO.h"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std::string_literals;
using namespace Scine;
using namespace molassembler;
using namespace DistanceGeometry;

/* TODO
 * - Prove EigenRefinementProblem is correct by verifying it against DlibRefinement
 * - Then prove EigenSIMDRefinementProblem is correct by verifying it against
 *   EigenRefinementProblem
 * - Then remove DlibRefinement and optimize either variant for speed,
 *   verifying them against each other
 *
 * - Planned tests:
 *   - Rotational and translational invariance tests for all Eigen-based
 *     refinement problems
 *   - Check of optimized vs reference implementations (this will be a
 *     cross-test between Eigen-based refinement problems)
 */

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

struct RefinementBaseData {
  DistanceBoundsMatrix distanceBounds;
  Eigen::MatrixXd embeddedPositions;
  std::vector<ChiralityConstraint> chiralConstraints;
  std::vector<DihedralConstraint> dihedralConstraints;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Eigen::MatrixXd squaredBounds() const {
    return distanceBounds.access().cwiseProduct(distanceBounds.access());
  }

  Eigen::VectorXd linearizeEmbeddedPositions() {
    return Eigen::Map<Eigen::VectorXd>(
      embeddedPositions.data(),
      embeddedPositions.cols() * embeddedPositions.rows()
    );
  }

  RefinementBaseData() = default;

  RefinementBaseData(const std::string& filename) {
    Molecule molecule = IO::read(filename);

    auto DGInfo = gatherDGInformation(molecule, DistanceGeometry::Configuration {});

    distanceBounds = DistanceBoundsMatrix {
      molecule,
      DGInfo.bounds
    };

    chiralConstraints = std::move(DGInfo.chiralityConstraints);
    dihedralConstraints = std::move(DGInfo.dihedralConstraints);

    auto distancesResult = distanceBounds.makeDistanceMatrix(randomnessEngine());
    if(!distancesResult) {
      BOOST_FAIL(distancesResult.error().message());
    }
    auto distances = distancesResult.value();

    MetricMatrix metric(distances);

    embeddedPositions = metric.embed();
  }
};

BOOST_AUTO_TEST_CASE( numericDifferentiationApproximatesGradients ) {
  using Vector = dlib::matrix<double, 0, 1>;

  Log::level = Log::Level::None;
  Log::particulars = {};

  // Generate several RefinementProblems and check their gradients
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("ez_stereocenters")
  ) {
    RefinementBaseData baseData {currentFilePath.string()};

    // Transfer to dlib
    const Vector dlibPositions = dlib::mat(
      baseData.linearizeEmbeddedPositions()
    );

    const dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      baseData.squaredBounds()
    );

    ErrorFunctionValue valueFunctor {
      squaredBounds,
      baseData.chiralConstraints,
      baseData.dihedralConstraints
    };

    ErrorFunctionGradient gradientFunctor {
      squaredBounds,
      baseData.chiralConstraints,
      baseData.dihedralConstraints
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

    valueFunctor.compressFourthDimension = true;
    gradientFunctor.compressFourthDimension = true;

    Vector compressedFiniteDifferenceGradient = dlib::derivative(valueFunctor)(dlibPositions);
    Vector compressedGradient = gradientFunctor(dlibPositions);

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
    RefinementBaseData baseData {currentFilePath.string()};

    // Transfer to dlib and vectorize
    ErrorFunctionValue::Vector referencePositions = dlib::mat(
      baseData.linearizeEmbeddedPositions()
    );

    dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      baseData.squaredBounds()
    );

    ErrorFunctionValue valueFunctor {
      squaredBounds,
      baseData.chiralConstraints,
      baseData.dihedralConstraints
    };
    valueFunctor.compressFourthDimension = true;

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
      Eigen::MatrixXd embeddedCopy = baseData.embeddedPositions;

      // Transform the embedded coordinates by rotation and translation
      for(unsigned i = 0; i < embeddedCopy.cols(); i++) {
        embeddedCopy.template block<3, 1>(0, i) = (
          rotationMatrix * embeddedCopy.template block<3, 1>(0, i)
        );
        embeddedCopy.template block<3, 1>(0, i) += translation;
      }

      // Transfer to dlib and vectorize
      ErrorFunctionValue::Vector vectorizedPositions = dlib::mat(
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
    RefinementBaseData baseData {currentFilePath.string()};

    // Vectorize the positions for use with dlib
    ErrorFunctionValue::Vector referencePositions = dlib::mat(
      baseData.linearizeEmbeddedPositions()
    );

    dlib::matrix<double, 0, 0> squaredBounds = dlib::mat(
      baseData.squaredBounds()
    );

    ErrorFunctionGradient gradientFunctor {
      squaredBounds,
      baseData.chiralConstraints,
      baseData.dihedralConstraints
    };

    gradientFunctor.compressFourthDimension = true;

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
      Eigen::MatrixXd rotatedRectangularPositions = baseData.embeddedPositions;
      Eigen::MatrixXd translatedRectangularPositions = baseData.embeddedPositions;

      // Transform the embedded coordinates by rotation and translation
      for(unsigned i = 0; i < rotatedRectangularPositions.cols(); i++) {
        rotatedRectangularPositions.template block<3, 1>(0, i) = (
          rotationMatrix * rotatedRectangularPositions.template block<3, 1>(0, i)
        );
        translatedRectangularPositions.template block<3, 1>(0, i) += translation;
      }

      // Vectorize the positions
      ErrorFunctionValue::Vector rotatedPositions = dlib::mat(
        static_cast<Eigen::VectorXd>(
          Eigen::Map<Eigen::VectorXd>(
            rotatedRectangularPositions.data(),
            rotatedRectangularPositions.cols() * rotatedRectangularPositions.rows()
          )
        )
      );

      ErrorFunctionValue::Vector translatedPositions = dlib::mat(
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

namespace traits {

template<typename FloatType>
struct Floating {};

template<>
struct Floating<double> {
  constexpr static double tolerance = 1e-10;
  constexpr static const char* name = "double";
};

template<>
struct Floating<float> {
  constexpr static double tolerance = 1e-5;
  constexpr static const char* name = "float";
};

constexpr double Floating<double>::tolerance;
constexpr const char* Floating<double>::name;
constexpr double Floating<float>::tolerance;
constexpr const char* Floating<float>::name;

} // namespace traits

Eigen::MatrixXd rotateAndTranslate(
  Eigen::MatrixXd positionMatrix,
  const Eigen::Matrix3d& rotationMatrix,
  const Eigen::Vector3d& translationVector
) {
  assert(positionMatrix.rows() == 4 || positionMatrix.rows() == 3);

  /* Transform the copy of embedded coordinates by rotation and translation
   * (only x,y,z coordinates, even if a fourth spatial coordinate is present)
   */
  for(unsigned i = 0; i < positionMatrix.cols(); ++i) {
    positionMatrix.template block<3, 1>(0, i) = (
      rotationMatrix * positionMatrix.template block<3, 1>(0, i)
    );
    positionMatrix.template block<3, 1>(0, i) += translationVector;
  }

  return positionMatrix;
}

template<typename EigenRefinementType>
Eigen::VectorXd rotateAndTranslateLinear(
  Eigen::VectorXd positionVector,
  const Eigen::Matrix3d& rotationMatrix,
  const Eigen::Vector3d& translationVector
) {
  constexpr unsigned dimensionality = EigenRefinementType::dimensions;

  assert(positionVector.size() % dimensionality == 0);
  const unsigned N = positionVector.size() / dimensionality;
  for(unsigned i = 0; i < N; ++i) {
    positionVector.segment<3>(dimensionality * i) = (
      rotationMatrix * positionVector.segment<3>(dimensionality * i)
      + translationVector
    );
  }

  return positionVector;
}

Eigen::VectorXd flatten(Eigen::MatrixXd matr) {
  return Eigen::Map<Eigen::VectorXd>(
    matr.data(),
    matr.cols() * matr.rows()
  );
}

template<typename EigenRefinementType>
struct RotationalTranslationalInvarianceTest {
  using PositionType = typename EigenRefinementType::VectorType;

  struct ErrorAndGradientComponent {
    using FunctionType = std::function<void(const PositionType&, double&, Eigen::VectorXd&)>;

    std::string name;
    FunctionType function;

    double referenceValue;
    Eigen::VectorXd referenceGradient;
    bool errorRotationalInvariance = true;
    bool errorTranslationalInvariance = true;
    bool gradientRotationalInvariance = true;
    bool gradientTranslationalInvariance = true;

    // Initialize name, function and reference value and gradient
    void initialize(std::string passName, FunctionType passFunction, const PositionType& referencePositions) {
      name = passName;
      function = passFunction;

      referenceValue = 0;
      referenceGradient.resize(referencePositions.size());
      referenceGradient.setZero();

      function(referencePositions, referenceValue, referenceGradient);
    }

    /* Propagate errorInvariance and gradientInvariance with new rotationally
     * and translationally transformed coordinates. The function is not called
     * again if invariance has been disproved for both error value and gradient.
     */
    void propagate(
      const Eigen::MatrixXd& positionMatrix,
      const Eigen::Matrix3d& rotationMatrix,
      const Eigen::Vector3d& translationVector
    ) {
      if(errorTranslationalInvariance || gradientTranslationalInvariance) {
        double translatedValue = 0;
        Eigen::VectorXd translatedGradient;
        translatedGradient.resize(positionMatrix.cols() * positionMatrix.rows());
        translatedGradient.setZero();

        PositionType translatedPositions = EigenRefinementType::conditionalDowncast(
          flatten(
            rotateAndTranslate(
              positionMatrix,
              Eigen::Matrix3d::Identity(),
              translationVector
            )
          )
        );

        function(translatedPositions, translatedValue, translatedGradient);

        if(
          errorTranslationalInvariance
          && (
            std::fabs(translatedValue - referenceValue)
            > traits::Floating<typename EigenRefinementType::FloatingPointType>::tolerance
          )
        ) {
          std::cout << "Error is not translationally invariant: translated = "
            << translatedValue << ", reference = " << referenceValue
            << ", difference: " << std::fabs(translatedValue - referenceValue)
            << "\n";
          errorTranslationalInvariance = false;
        }

        if(gradientTranslationalInvariance) {
          if(
            !translatedGradient.isApprox(
              referenceGradient,
              traits::Floating<typename EigenRefinementType::FloatingPointType>::tolerance
            )
          ) {
            std::cout << "Translated gradient: " << translatedGradient.transpose() << "\n";
            std::cout << "Reference gradient : " << referenceGradient.transpose() << "\n";
            std::cout << "Trans. - reference : " << (translatedGradient - referenceGradient).transpose() << "\n";
            gradientTranslationalInvariance = false;
          }
        }
      }

      if(errorRotationalInvariance || gradientRotationalInvariance) {
        double rotatedValue = 0;
        Eigen::VectorXd rotatedGradient;
        rotatedGradient.resize(positionMatrix.cols() * positionMatrix.rows());
        rotatedGradient.setZero();

        PositionType rotatedPositions = EigenRefinementType::conditionalDowncast(
          flatten(
            rotateAndTranslate(
              positionMatrix,
              rotationMatrix,
              Eigen::Vector3d::Zero()
            )
          )
        );

        function(rotatedPositions, rotatedValue, rotatedGradient);

        if(
          errorRotationalInvariance
          && (
            std::fabs(rotatedValue - referenceValue)
            > traits::Floating<typename EigenRefinementType::FloatingPointType>::tolerance
          )
        ) {
          std::cout << "Error is not rotationally invariant: rotated = "
            << rotatedValue << ", reference = " << referenceValue
            << ", difference: " << std::fabs(rotatedValue - referenceValue)
            << "\n";
          errorRotationalInvariance = false;
        }

        if(gradientRotationalInvariance) {
          // Rotate reference accordingly
          auto rotatedReferenceGradient = rotateAndTranslateLinear<EigenRefinementType>(
            referenceGradient,
            rotationMatrix,
            Eigen::Vector3d::Zero()
          );

          if(
            !rotatedGradient.isApprox(
              rotatedReferenceGradient,
              traits::Floating<typename EigenRefinementType::FloatingPointType>::tolerance
            )
          ) {
            std::cout << "Rotated   gradient: " << rotatedGradient.transpose() << "\n";
            std::cout << "Rot. ref. gradient: " << rotatedReferenceGradient.transpose() << "\n";
            std::cout << "Rot. - reference  : " << (rotatedGradient - rotatedReferenceGradient).transpose() << "\n";
            gradientRotationalInvariance = false;
          }
        }
      }
    }
  };

  static bool value() {
    bool pass = true;

    for(
      const boost::filesystem::path& currentFilePath :
      boost::filesystem::recursive_directory_iterator("ez_stereocenters")
    ) {
      RefinementBaseData baseData {currentFilePath.string()};

      EigenRefinementType refinementFunctor {
        baseData.squaredBounds(),
        baseData.chiralConstraints,
        baseData.dihedralConstraints
      };

      using PositionsType = typename EigenRefinementType::VectorType;
      using FloatingPointType = typename EigenRefinementType::FloatingPointType;

      const PositionsType referencePositions = baseData.linearizeEmbeddedPositions().template cast<FloatingPointType>();

      std::array<ErrorAndGradientComponent, 4> components;

      auto makeFunctor = [&refinementFunctor](auto memFn) {
        return [&refinementFunctor, memFn](
          const PositionsType& positions,
          double& value,
          Eigen::VectorXd& gradient
        ) { memFn(refinementFunctor, positions, value, gradient); };
      };

      components[0].initialize(
        "distance",
        makeFunctor(std::mem_fn(&EigenRefinementType::distanceContributions)),
        referencePositions
      );

      components[1].initialize(
        "chiral",
        makeFunctor(std::mem_fn(&EigenRefinementType::chiralContributions)),
        referencePositions
      );

      components[2].initialize(
        "dihedral",
        makeFunctor(std::mem_fn(&EigenRefinementType::dihedralContributions)),
        referencePositions
      );

      components[3].initialize(
        "fourth dimension",
        makeFunctor(std::mem_fn(&EigenRefinementType::fourthDimensionContributions)),
        referencePositions
      );

      for(unsigned testNum = 0; testNum < 100; ++testNum) {
        const Eigen::Matrix3d rotationMatrix = Eigen::Quaterniond::UnitRandom().toRotationMatrix();
        const Eigen::Vector3d translationVector = Eigen::Vector3d::Random();

        for(auto& component : components) {
          component.propagate(baseData.embeddedPositions, rotationMatrix, translationVector);
        }
      }

      bool filePasses = temple::all_of(
        components,
        [](const ErrorAndGradientComponent& component) -> bool {
          return (
            component.errorTranslationalInvariance
            && component.errorRotationalInvariance
            && component.gradientTranslationalInvariance
            && component.gradientRotationalInvariance
          );
        }
      );

      if(!filePasses) {
        std::cout << "For file " << currentFilePath.string() << ":\n";
        for(const auto& component : components) {
          std::string compound = "- " + component.name + ": error invariance (rot: ";

          auto addBooleanStr = [&compound](bool value) {
            if(value) {
              compound += "yes";
            } else {
              compound += "NO";
            }
          };

          addBooleanStr(component.errorRotationalInvariance);
          compound +=", trans: ";
          addBooleanStr(component.errorTranslationalInvariance);
          compound += "), gradient invariance (rot: ";
          addBooleanStr(component.gradientRotationalInvariance);
          compound +=", trans: ";
          addBooleanStr(component.gradientTranslationalInvariance);

          compound += ")\n";

          std::cout << compound;
        }
      }

      pass &= filePasses;
    }

    if(!pass) {
      std::cout << "\n\n";
    }

    return pass;
  }
};

BOOST_AUTO_TEST_CASE(RefinementProblemRotationalTranslationalInvariance) {
  using EigenRefinementTypeVariations = std::tuple<
    EigenRefinementProblem<4, double>,
    EigenRefinementProblem<4, float>,
    EigenSIMDRefinementProblem<4, double>,
    EigenSIMDRefinementProblem<4, float>
  >;

  constexpr unsigned variations = std::tuple_size<EigenRefinementTypeVariations>::value;

  std::array<std::string, variations> eigenRefinementNames {
    "EigenRefinementProblem<4, double>",
    "EigenRefinementProblem<4, float>",
    "EigenSIMDRefinementProblem<4, double>",
    "EigenSIMDRefinementProblem<4, float>"
  };

  auto passes = temple::TupleType::map<EigenRefinementTypeVariations, RotationalTranslationalInvarianceTest>();

  for(unsigned i = 0; i < variations; ++i) {
    BOOST_CHECK_MESSAGE(
      passes.at(i),
      "Rotational and translational invariance fail for " << eigenRefinementNames.at(i)
    );
  }
}
