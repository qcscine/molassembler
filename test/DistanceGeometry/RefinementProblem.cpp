/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "Eigen/Geometry"
#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Functor.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/constexpr/FloatingPointComparison.h"
#include "Molassembler/Temple/constexpr/Numeric.h"
#include "Molassembler/Temple/constexpr/TupleType.h"
#include "Molassembler/Temple/constexpr/TupleTypePairs.h"

#include "Molassembler/Modeling/CommonTrig.h"
#include "Molassembler/DistanceGeometry/ConformerGeneration.h"
#include "Molassembler/DistanceGeometry/MetricMatrix.h"
#include "Molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "Molassembler/DistanceGeometry/EigenRefinement.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/IO.h"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std::string_literals;
using namespace Scine::Molassembler;
using namespace DistanceGeometry;

struct RefinementBaseData {
  DistanceBoundsMatrix distanceBounds;
  Eigen::MatrixXd embeddedPositions;
  std::vector<ChiralConstraint> chiralConstraints;
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

    auto DgInfo = gatherDGInformation(molecule, DistanceGeometry::Configuration {});

    distanceBounds = DistanceBoundsMatrix {
      molecule.graph().inner(),
      DgInfo.bounds
    };

    chiralConstraints = std::move(DgInfo.chiralConstraints);
    dihedralConstraints = std::move(DgInfo.dihedralConstraints);

    auto distancesResult = distanceBounds.makeDistanceMatrix(randomnessEngine());
    if(!distancesResult) {
      BOOST_FAIL(distancesResult.error().message());
    }
    auto distances = distancesResult.value();

    MetricMatrix metric(distances);

    embeddedPositions = metric.embed();
  }
};

namespace Traits {

template<typename FloatType>
struct Floating {};

template<>
struct Floating<double> {
  constexpr static double tolerance = 1e-8;
};

template<>
struct Floating<float> {
  constexpr static double tolerance = 1e-4;
};

constexpr double Floating<double>::tolerance;
constexpr double Floating<float>::tolerance;

} // namespace Traits

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
auto rotateAndTranslateLine(
  const typename EigenRefinementType::VectorType& positionVector,
  const Eigen::Matrix3d& rotationMatrix,
  const Eigen::Vector3d& translationVector
) -> typename EigenRefinementType::VectorType {
  using FloatType = typename RefinementTraits<EigenRefinementType>::FloatingPointType;
  using DimensionalityConstant = typename RefinementTraits<EigenRefinementType>::DimensionalityConstant;
  constexpr unsigned dimensionality = DimensionalityConstant::value;

  typename EigenRefinementType::VectorType transformedPositions = positionVector;

  assert(positionVector.size() % dimensionality == 0);
  const unsigned N = positionVector.size() / dimensionality;
  for(unsigned i = 0; i < N; ++i) {
    Eigen::Vector3d position = positionVector.template segment<3>(dimensionality * i).template cast<double>();

    transformedPositions.template segment<3>(dimensionality * i) = (
      rotationMatrix * position + translationVector
    ).template cast<FloatType>();
  }

  return transformedPositions;
}

Eigen::VectorXd flatten(Eigen::MatrixXd matr) {
  return Eigen::Map<Eigen::VectorXd>(
    matr.data(),
    matr.cols() * matr.rows()
  );
}

template<typename EigenRefinementType>
struct RotationalTranslationalInvarianceTest {
  using FloatType = typename EigenRefinementType::FloatingPointType;
  using PositionType = typename EigenRefinementType::VectorType;

  struct ErrorAndGradientComponent {
    using FunctionType = std::function<void(const PositionType&, FloatType&, PositionType&)>;

    std::string name;
    FunctionType function;

    FloatType referenceValue;
    PositionType referenceGradient;
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
        FloatType translatedValue = 0;
        PositionType translatedGradient;
        translatedGradient.resize(positionMatrix.cols() * positionMatrix.rows());
        translatedGradient.setZero();

        PositionType translatedPositions = flatten(
          rotateAndTranslate(
            positionMatrix,
            Eigen::Matrix3d::Identity(),
            translationVector
          )
        ).cast<FloatType>();

        function(translatedPositions, translatedValue, translatedGradient);

        const double translatedDifference = std::fabs(
          static_cast<double>(translatedValue)
          - static_cast<double>(referenceValue)
        );

        if(
          errorTranslationalInvariance
          && translatedDifference > Traits::Floating<FloatType>::tolerance
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
              Traits::Floating<FloatType>::tolerance
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
        FloatType rotatedValue = 0;
        PositionType rotatedGradient;
        rotatedGradient.resize(positionMatrix.cols() * positionMatrix.rows());
        rotatedGradient.setZero();

        PositionType rotatedPositions = flatten(
          rotateAndTranslate(
            positionMatrix,
            rotationMatrix,
            Eigen::Vector3d::Zero()
          )
        ).cast<FloatType>();

        function(rotatedPositions, rotatedValue, rotatedGradient);

        const double rotatedDifference = std::fabs(
          static_cast<double>(rotatedValue)
          - static_cast<double>(referenceValue)
        );

        if(
          errorRotationalInvariance
          && rotatedDifference > Traits::Floating<FloatType>::tolerance

        ) {
          std::cout << "Error is not rotationally invariant: rotated = "
            << rotatedValue << ", reference = " << referenceValue
            << ", difference: " << rotatedDifference
            << "\n";
          errorRotationalInvariance = false;
        }

        if(gradientRotationalInvariance) {
          // Rotate reference accordingly
          auto rotatedReferenceGradient = rotateAndTranslateLine<EigenRefinementType>(
            referenceGradient,
            rotationMatrix,
            Eigen::Vector3d::Zero()
          );

          if(
            !rotatedGradient.isApprox(
              rotatedReferenceGradient,
              Traits::Floating<FloatType>::tolerance
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
          FloatingPointType& value,
          PositionsType& gradient
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

      bool filePasses = Temple::all_of(
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
    EigenRefinementProblem<4, double, false>,
    //EigenRefinementProblem<4, float, false>,
    EigenRefinementProblem<4, double, true>
    //EigenRefinementProblem<4, float, true>
  >;

  constexpr unsigned variations = std::tuple_size<EigenRefinementTypeVariations>::value;

  std::array<std::string, variations> eigenRefinementNames {
    "EigenRefinementProblem<4, double, false>",
    //"EigenRefinementProblem<4, float, false>",
    "EigenRefinementProblem<4, double, true>"
    // "EigenRefinementProblem<4, float, true>"
  };

  auto passes = Temple::Tuples::map<EigenRefinementTypeVariations, RotationalTranslationalInvarianceTest>();

  for(unsigned i = 0; i < variations; ++i) {
    BOOST_CHECK_MESSAGE(
      passes.at(i),
      "Rotational and translational invariance fail for " << eigenRefinementNames.at(i)
    );
  }
}

template<typename RefinementT, typename RefinementU>
struct CompareImplementations {
  using TraitsT = RefinementTraits<RefinementT>;
  using TraitsU = RefinementTraits<RefinementU>;

  using PositionsT = typename RefinementT::VectorType;
  using FloatingPointT = typename TraitsT::FloatingPointType;
  using PositionsU = typename RefinementU::VectorType;
  using FloatingPointU = typename TraitsU::FloatingPointType;

  using SmallerFPType = std::conditional_t<
    sizeof(FloatingPointT) < sizeof(FloatingPointU),
    FloatingPointT,
    FloatingPointU
  >;

  using LargerFPType = std::conditional_t<
    sizeof(FloatingPointT) < sizeof(FloatingPointU),
    FloatingPointU,
    FloatingPointT
  >;

  static bool value() {
    return Temple::all_of(
      boost::filesystem::recursive_directory_iterator("ez_stereocenters"),
      [](const boost::filesystem::path& currentFilePath) -> bool {
        RefinementBaseData baseData {currentFilePath.string()};

        RefinementT tFunctor {
          baseData.squaredBounds(),
          baseData.chiralConstraints,
          baseData.dihedralConstraints
        };

        RefinementU uFunctor {
          baseData.squaredBounds(),
          baseData.chiralConstraints,
          baseData.dihedralConstraints
        };

        PositionsT tPositions = baseData.linearizeEmbeddedPositions().template cast<FloatingPointT>();
        PositionsU uPositions = baseData.linearizeEmbeddedPositions().template cast<FloatingPointU>();

        using FunctionTypeT = std::function<void(const PositionsT&, FloatingPointT&, PositionsT&)>;
        using FunctionTypeU = std::function<void(const PositionsU&, FloatingPointU&, PositionsU&)>;

        std::vector<
          std::tuple<std::string, FunctionTypeT, FunctionTypeU>
        > components {
          {
            "distance",
            [&tFunctor](const PositionsT& positions, FloatingPointT& value, PositionsT& gradients) {
              tFunctor.distanceContributions(positions, value, gradients);
            },
            [&uFunctor](const PositionsU& positions, FloatingPointU& value, PositionsU& gradients) {
              uFunctor.distanceContributions(positions, value, gradients);
            }
          },
          {
            "chiral",
            [&tFunctor](const PositionsT& positions, FloatingPointT& value, PositionsT& gradients) {
              tFunctor.chiralContributions(positions, value, gradients);
            },
            [&uFunctor](const PositionsU& positions, FloatingPointU& value, PositionsU& gradients) {
              uFunctor.chiralContributions(positions, value, gradients);
            }
          },
          {
            "dihedral",
            [&tFunctor](const PositionsT& positions, FloatingPointT& value, PositionsT& gradients) {
              tFunctor.dihedralContributions(positions, value, gradients);
            },
            [&uFunctor](const PositionsU& positions, FloatingPointU& value, PositionsU& gradients) {
              uFunctor.dihedralContributions(positions, value, gradients);
            }
          },
          {
            "fourth dimension",
            [&tFunctor](const PositionsT& positions, FloatingPointT& value, PositionsT& gradients) {
              tFunctor.fourthDimensionContributions(positions, value, gradients);
            },
            [&uFunctor](const PositionsU& positions, FloatingPointU& value, PositionsU& gradients) {
              uFunctor.fourthDimensionContributions(positions, value, gradients);
            }
          },
        };

        auto passesComparisonMap = Temple::map(
          components,
          [&](const auto& comparisonTuple) -> bool {
            FloatingPointT tError = 0;
            FloatingPointU uError = 0;
            PositionsT tGradients;
            PositionsU uGradients;
            tGradients.resize(tPositions.size());
            tGradients.setZero();
            uGradients.resize(uPositions.size());
            uGradients.setZero();

            std::get<1>(comparisonTuple)(tPositions, tError, tGradients);
            std::get<2>(comparisonTuple)(uPositions, uError, uGradients);

            bool pass = true;

            const double errorAbsDifference = std::fabs(
              static_cast<double>(tError)
              - static_cast<double>(uError)
            );

            if(errorAbsDifference > Traits::Floating<SmallerFPType>::tolerance) {
              std::cout << "Error values for component "
                << std::get<0>(comparisonTuple) << " do not match between "
                << RefinementT::name() << " and " << RefinementU::name() << " for "
                << currentFilePath.string() << ": "
                << tError << " != " << uError << ", difference = "
                << errorAbsDifference << "\n";

              pass = false;
            }

            if(
              !tGradients.template cast<LargerFPType>().isApprox(
                uGradients.template cast<LargerFPType>(),
                Traits::Floating<SmallerFPType>::tolerance
              )
            ) {
              std::cout << "Gradients for component "
                << std::get<0>(comparisonTuple) << " do not match between "
                << RefinementT::name() << " and " << RefinementU::name() << "  for "
                << currentFilePath.string() << ":\n  "
                << tGradients.transpose() << "\n  "
                << uGradients.transpose() << "\n";

              pass = false;
            }

            return pass;
          }
        );

        return Temple::all_of(passesComparisonMap);
      }
    );
  }
};

BOOST_AUTO_TEST_CASE(RefinementProblemEquivalence) {
  using DoubleRefinementTypeVariations = std::tuple<
    EigenRefinementProblem<4, double, false>,
    EigenRefinementProblem<4, double, true>
  >;

  auto doublePasses = Temple::Tuples::mapAllPairs<DoubleRefinementTypeVariations, CompareImplementations>();

  BOOST_CHECK_MESSAGE(
    Temple::all_of(doublePasses),
    "Not all refinement template argument of double variations match pair-wise!"
  );

  using FloatRefinementTypeVariations = std::tuple<
    EigenRefinementProblem<4, float, false>,
    EigenRefinementProblem<4, float, true>
  >;

  auto floatPasses = Temple::Tuples::mapAllPairs<FloatRefinementTypeVariations, CompareImplementations>();

  BOOST_CHECK_MESSAGE(
    Temple::all_of(floatPasses),
    "Not all refinement template argument of float variations match pair-wise!"
  );
}
