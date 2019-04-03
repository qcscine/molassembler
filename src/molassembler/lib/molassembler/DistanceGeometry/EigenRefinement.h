/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Eigen refinement implementation
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_EIGEN_REFINEMENT_PROBLEM_H
#define INCLUDE_MOLASSEMBLER_DG_EIGEN_REFINEMENT_PROBLEM_H

#include <Eigen/Core>
#include "molassembler/Types.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"

namespace Scine {
namespace molassembler {
namespace DistanceGeometry {

template<unsigned dimensionality, typename FloatType>
struct EigenRefinementProblem {
  static_assert(
    dimensionality == 3 || dimensionality == 4,
    "RefinementErrorFunction is only suitable for three or four dimensions"
  );

//!@name Public types
//!@{
  //! Vector layout of positions
  using VectorType = Eigen::Matrix<FloatType, 1, Eigen::Dynamic>;

  //! Three dimensional vector
  using ThreeDimensionalVector = Eigen::Matrix<FloatType, 3, 1>;
  //! Three or four dimensional vector (depending on dimensionality)
  using FullDimensionalVector = Eigen::Matrix<FloatType, dimensionality, 1>;

  //! Three-row dynamic column matrix
  using ThreeDimensionalMatrixType = Eigen::Matrix<FloatType, 3, Eigen::Dynamic>;
  //! Three or four-row dynamic column matrix
  using FullDimensionalMatrixType = Eigen::Matrix<FloatType, dimensionality, Eigen::Dynamic>;
  //! Template argument specifying floating-point type
  using FloatingPointType = FloatType;
//!@}

//!@name Static members
//!@{
  constexpr static unsigned dimensions = dimensionality;
//!@}

  template<
    typename EigenType,
    typename FPType = FloatType,
    typename std::enable_if_t<std::is_same<FPType, double>::value, int> = 0
  > static auto conditionalDowncast(EigenType eigenTypeValue) {
    return eigenTypeValue;
  }

  template<
    typename EigenType,
    typename FPType = FloatType,
    typename std::enable_if_t<!std::is_same<FPType, double>::value, int> = 0
  > static auto conditionalDowncast(EigenType eigenTypeValue) {
    return eigenTypeValue.template cast<FloatType>().eval();
  }

  template<
    typename EigenType,
    typename FPType = FloatType,
    typename std::enable_if_t<std::is_same<FPType, double>::value, int> = 0
  > static auto conditionalUpcast(EigenType eigenTypeValue) {
    return eigenTypeValue;
  }

  template<
    typename EigenType,
    typename FPType = FloatType,
    typename std::enable_if_t<!std::is_same<FPType, double>::value, int> = 0
  > static auto conditionalUpcast(EigenType eigenTypeValue) {
    return eigenTypeValue.template cast<double>().eval();
  }

//!@name Static member functions
//!@{
  //! Copy a full dimensional position of an atom
  inline static FullDimensionalVector getPosition(
    const VectorType& positions,
    const AtomIndex index
  ) {
    return positions.template segment<dimensionality>(dimensionality * index);
  }

  //! Copy a three-dimensional position of an atom
  inline static ThreeDimensionalVector getPosition3D(
    const VectorType& positions,
    const AtomIndex index
  ) {
    return positions.template segment<3>(dimensionality * index);
  }

  //! Calculate the average full-dimensional position of several atoms
  static FullDimensionalVector getAveragePosition(
    const VectorType& positions,
    const ChiralityConstraint::AtomListType& atomList
  ) {
    if(atomList.size() == 1) {
      return getPosition(positions, atomList.front());
    }

    FullDimensionalVector sum = FullDimensionalVector::Zero();
    for(const AtomIndex i : atomList) {
      sum += positions.template segment<dimensionality>(dimensionality * i);
    }

    return sum / atomList.size();
  }

  //! Calculate the average three-dimensional position of several atoms
  static ThreeDimensionalVector getAveragePosition3D(
    const VectorType& positions,
    const ChiralityConstraint::AtomListType& atomList
  ) {
    if(atomList.size() == 1) {
      return getPosition3D(positions, atomList.front());
    }

    ThreeDimensionalVector sum = ThreeDimensionalVector::Zero();
    for(const AtomIndex i : atomList) {
      sum += positions.template segment<3>(dimensionality * i);
    }

    return sum / atomList.size();
  }
//!@}

//!@name Public members
//!@{
  Eigen::MatrixXd squaredBounds;
  std::vector<ChiralityConstraint> chiralConstraints;
  std::vector<DihedralConstraint> dihedralConstraints;
  bool compressFourthDimension = false;
//!@}

//!@name Signaling members
//!@{
  mutable double proportionChiralConstraintsCorrectSign = 0.0;
//!@}

  explicit EigenRefinementProblem(
    Eigen::MatrixXd passSquaredBounds,
    std::vector<ChiralityConstraint> passChiralityConstraints,
    std::vector<DihedralConstraint> passDihedralConstraints
  ) : squaredBounds(std::move(passSquaredBounds)),
      chiralConstraints(std::move(passChiralityConstraints)),
      dihedralConstraints(std::move(passDihedralConstraints))
  {}

  /*!
   * @brief Adds distance error and gradient contributions
   */
  void distanceContributions(
    const VectorType& positions,
    double& error,
    Eigen::VectorXd& gradients
  ) const {
    assert(gradients.size() == positions.size());
    const unsigned N = positions.size() / dimensionality;

    for(unsigned i = 0; i < N - 1; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        const double lowerBoundSquared = squaredBounds(j, i);
        const double upperBoundSquared = squaredBounds(i, j);
        assert(lowerBoundSquared <= upperBoundSquared);

        // For both
        const FullDimensionalVector positionDifference = (
          positions.template segment<dimensionality>(dimensionality * i)
          - positions.template segment<dimensionality>(dimensionality * j)
        );

        const double diffLength = positionDifference.norm();

        // Upper term
        const double upperTerm = diffLength / upperBoundSquared - 1;

        assert(!( // Ensure early return logic is correct
          upperTerm > 0
          && ( // lowerTerm
            2 * lowerBoundSquared / (
              lowerBoundSquared + diffLength
            ) - 1
          ) > 0
        ));

        if(upperTerm > 0) {
          error += upperTerm * upperTerm;

          const auto f = conditionalUpcast(
            FullDimensionalVector {
              4 * positionDifference * upperTerm / upperBoundSquared
            }
          );

          gradients.template segment<dimensionality>(dimensionality * i) += f;
          gradients.template segment<dimensionality>(dimensionality * j) -= f;

          /* if the upper term is set, there is no way the lower term needs to be
           * applied -> early return
           */
          return;
        }

        // Lower term
        const double quotient = lowerBoundSquared + diffLength;
        const double lowerTerm = 2 * lowerBoundSquared / quotient - 1;

        if(lowerTerm > 0) {
          error += lowerTerm * lowerTerm;

          const auto g = conditionalUpcast(
            FullDimensionalVector {
              8 * lowerBoundSquared * positionDifference * lowerTerm / (
                quotient * quotient
              )
            }
          );

          /* We use -= because the lower term needs the position vector
           * difference (j - i), so we reuse positionDifference and just subtract
           * from the gradient instead of adding to it
           */
          gradients.template segment<dimensionality>(dimensionality * i) -= g;
          gradients.template segment<dimensionality>(dimensionality * j) += g;
        }
      }
    }
  }

  /*!
   * @brief Adds chiral error and gradient contributions
   */
  void chiralContributions(const VectorType& positions, double& error, Eigen::VectorXd& gradients) const {
    unsigned nonZeroChiralityConstraints = 0;
    unsigned incorrectNonZeroChiralityConstraints = 0;

    for(const auto& constraint : chiralConstraints) {
      const ThreeDimensionalVector alpha = getAveragePosition3D(positions, constraint.sites[0]);
      const ThreeDimensionalVector beta = getAveragePosition3D(positions, constraint.sites[1]);
      const ThreeDimensionalVector gamma = getAveragePosition3D(positions, constraint.sites[2]);
      const ThreeDimensionalVector delta = getAveragePosition3D(positions, constraint.sites[3]);

      const ThreeDimensionalVector alphaMinusDelta = alpha - delta;
      const ThreeDimensionalVector betaMinusDelta = beta - delta;
      const ThreeDimensionalVector gammaMinusDelta = gamma - delta;

      const double volume = (alphaMinusDelta).dot(
        betaMinusDelta.cross(gammaMinusDelta)
      );

      if(std::fabs(constraint.lower + constraint.upper) > 1e-4) {
        ++nonZeroChiralityConstraints;
        if(
          ( volume < 0 && constraint.lower > 0)
          || (volume > 0 && constraint.lower < 0)
        ) {
          ++incorrectNonZeroChiralityConstraints;
        }
      }

      // Upper bound term
      const double upperTerm = volume - constraint.upper;

      assert(!( // Ensure early continue logic is correct
        upperTerm > 0
        && constraint.lower - volume > 0
      ));

      if(upperTerm > 0) {
        error += upperTerm * upperTerm;
      }

      // Lower bound term
      const double lowerTerm = constraint.lower - volume;
      if(lowerTerm > 0) {
        error += lowerTerm * lowerTerm;
      }

      const double factor = 2 * (
        std::max(0.0, upperTerm)
        - std::max(0.0, lowerTerm)
      );

      /* Make sure computing all cross products is worth it by checking that
       * one of both terms is actually greater than 0
       */
      if(factor != 0.0) {
        // Specific to deltaI only but still repeated
        const ThreeDimensionalVector alphaMinusGamma = alpha - gamma;

        const ThreeDimensionalVector betaMinusGamma = beta - gamma;

        const auto iContribution = conditionalUpcast(
          ThreeDimensionalVector {
            factor * betaMinusDelta.cross(gammaMinusDelta)
            / constraint.sites[0].size()
          }
        );

        const auto jContribution = conditionalUpcast(
          ThreeDimensionalVector {
            factor * gammaMinusDelta.cross(alphaMinusDelta)
            / constraint.sites[1].size()
          }
        );

        const auto kContribution = conditionalUpcast(
          ThreeDimensionalVector {
            factor * alphaMinusDelta.cross(betaMinusDelta)
            / constraint.sites[2].size()
          }
        );

        const auto lContribution = conditionalUpcast(
          ThreeDimensionalVector {
            factor * betaMinusGamma.cross(alphaMinusGamma)
            / constraint.sites[3].size()
          }
        );

        for(const auto alphaI : constraint.sites[0]) {
          gradients.template segment<3>(dimensionality * alphaI) += iContribution;
        }

        for(const auto betaI : constraint.sites[1]) {
          gradients.template segment<3>(dimensionality * betaI) += jContribution;
        }

        for(const auto gammaI : constraint.sites[2]) {
          gradients.template segment<3>(dimensionality * gammaI) += kContribution;
        }

        for(const auto deltaI : constraint.sites[3]) {
          gradients.template segment<3>(dimensionality * deltaI) += lContribution;
        }
      }
    }

    // Set signaling member
    if(nonZeroChiralityConstraints == 0) {
      proportionChiralConstraintsCorrectSign = 1;
    } else {
      proportionChiralConstraintsCorrectSign = static_cast<double>(
        nonZeroChiralityConstraints - incorrectNonZeroChiralityConstraints
      ) / nonZeroChiralityConstraints;
    }
  }

  /*!
   * @brief Adds dihedral error and gradient contributions
   */
  void dihedralContributions(const VectorType& positions, double& error, Eigen::VectorXd& gradient) const {
    assert(gradient.size() == positions.size());
    constexpr double reductionFactor = 1.0 / 10;

    for(const DihedralConstraint& constraint : dihedralConstraints) {
      const ThreeDimensionalVector alpha = getAveragePosition3D(positions, constraint.sites[0]);
      const ThreeDimensionalVector beta = getAveragePosition3D(positions, constraint.sites[1]);
      const ThreeDimensionalVector gamma = getAveragePosition3D(positions, constraint.sites[2]);
      const ThreeDimensionalVector delta = getAveragePosition3D(positions, constraint.sites[3]);

      const ThreeDimensionalVector f = alpha - beta;
      const ThreeDimensionalVector g = beta - gamma;
      const ThreeDimensionalVector h = delta - gamma;

      ThreeDimensionalVector a = f.cross(g);
      ThreeDimensionalVector b = h.cross(g);

      // All of the permutations yield identical expressions within atan2 aside from -g
      double phi = std::atan2(
        a.cross(b).dot(
          -g.normalized()
        ),
        a.dot(b)
      );

      const double constraintSumHalved = (constraint.lower + constraint.upper) / 2;

      if(phi < constraintSumHalved - M_PI) {
        phi += 2 * M_PI;
      } else if(phi > constraintSumHalved + M_PI) {
        phi -= 2 * M_PI;
      }

      const double w_phi = phi - constraintSumHalved;

      double h_phi = std::fabs(w_phi) - (constraint.upper - constraint.lower) / 2;

      /* "Apply the max function": If h <= 0, then the max function yields zero,
       * so we can skip this constraint
       */
      if(h_phi <= 0) {
        continue;
      }

      // Error contribution
      error += h_phi * h_phi * reductionFactor;

      // Multiply with sgn (w)
      h_phi *= static_cast<int>(0.0 < w_phi) - static_cast<int>(w_phi < 0.0);

      // Multiply in the factors 2 * 1/10 -> 1/5
      h_phi *= 2 * reductionFactor;

      // Precompute some reused expressions
      const double gLength = g.norm();
      if(gLength == 0) {
        throw std::runtime_error("Encountered zero-length g vector in dihedral errors");
      }

      const double aLengthSq = a.squaredNorm();
      if(aLengthSq == 0) {
        throw std::runtime_error("Encountered zero-length a vector in dihedral gradient contributions");
      }
      a /= aLengthSq;
      const double bLengthSq = b.squaredNorm();
      if(bLengthSq == 0) {
        throw std::runtime_error("Encountered zero-length b vector in dihedral gradient contributions");
      }
      b /= bLengthSq;
      const double fDotG = f.dot(g);
      const double gDotH = g.dot(h);

      const auto iContribution = conditionalUpcast(
        ThreeDimensionalVector {
          -(h_phi / constraint.sites[0].size()) * gLength * a
        }
      );

      const auto jContribution = conditionalUpcast(
        ThreeDimensionalVector {
          (h_phi / constraint.sites[1].size()) * (
            (gLength + fDotG / gLength) * a
            - (gDotH / gLength) * b
          )
        }
      );

      const auto kContribution = conditionalUpcast(
        ThreeDimensionalVector {
          (h_phi / constraint.sites[2].size()) * (
            (gDotH / gLength - gLength) * b
            - (fDotG / gLength) * a
          )
        }
      );

      const auto lContribution = conditionalUpcast(
        ThreeDimensionalVector {
          (h_phi / constraint.sites[3].size()) * gLength * b
        }
      );

      if(
        !iContribution.array().isFinite().all()
        || !jContribution.array().isFinite().all()
        || !kContribution.array().isFinite().all()
        || !lContribution.array().isFinite().all()
      ) {
        throw std::out_of_range("Encountered non-finite dihedral contributions");
      }

      for(const AtomIndex alphaConstitutingIndex : constraint.sites[0]) {
        gradient.template segment<3>(dimensionality * alphaConstitutingIndex) += iContribution;
      }

      for(const AtomIndex betaConstitutingIndex: constraint.sites[1]) {
        gradient.template segment<3>(dimensionality * betaConstitutingIndex) += jContribution;
      }

      for(const AtomIndex gammaConstitutingIndex : constraint.sites[2]) {
        gradient.template segment<3>(dimensionality * gammaConstitutingIndex) += kContribution;
      }

      for(const AtomIndex deltaConstitutingIndex : constraint.sites[3]) {
        gradient.template segment<3>(dimensionality * deltaConstitutingIndex) += lContribution;
      }
    }
  }

  /*!
   * @brief Adds fourth dimension error and gradient contributions
   */
  void fourthDimensionContributions(const VectorType& positions, double& error, Eigen::VectorXd& gradient) const {
    assert(gradient.size() == positions.size());
    // Collect all fourth dimension values, square and sum, add to error
    if(dimensionality == 4 && compressFourthDimension) {
      const unsigned N = positions.size() / dimensionality;
      for(unsigned i = 0; i < N; ++i) {
        const FloatType fourthDimensionalValue = positions(4 * i + 3);
        error += static_cast<double>(fourthDimensionalValue * fourthDimensionalValue);
        gradient(4 * i + 3) += static_cast<double>(2 * fourthDimensionalValue);
      }
    }
  }

  /*!
   * @brief Calculates the error value and gradient for all contributions
   * @param parameters The linearized positions of all particles
   * @pre @p parameters must be evenly divisible by the dimensionality
   * @post Error is stored in @p value and gradient in @p gradient
   */
  void operator() (const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradient) const {
    assert(gradient.size() == parameters.size());
    assert(parameters.size() % dimensionality == 0);

    value = 0;
    gradient.setZero();

    /* Own downcasted parameters if FloatType is float, reference them
     * otherwise
     */
    std::conditional_t<
      std::is_same<FloatType, float>::value,
      VectorType,
      const VectorType&
    > possiblyDowncastedParameters = parameters;

    distanceContributions(possiblyDowncastedParameters, value, gradient);
    chiralContributions(possiblyDowncastedParameters, value, gradient);
    dihedralContributions(possiblyDowncastedParameters, value, gradient);
    fourthDimensionContributions(possiblyDowncastedParameters, value, gradient);
  }

  double calculateProportionChiralConstraintsCorrectSign(const Eigen::VectorXd& positions) const {
    unsigned nonZeroChiralityConstraints = 0;
    unsigned incorrectNonZeroChiralityConstraints = 0;

    for(const auto& constraint : chiralConstraints) {
      /* Make sure that chirality constraints meant to cause coplanarity of
       * four indices aren't counted here - Their volume bounds are
       * usually so tight that they might be slightly outside in any given
       * refinement step. If they are outside their target values by a large
       * amount, that is caught by finalStructureAcceptable anyway.
       */
      if(std::fabs(constraint.lower + constraint.upper) > 1e-4) {
        nonZeroChiralityConstraints += 1;

        const ThreeDimensionalVector delta = getAveragePosition3D(positions, constraint.sites[3]);

        const double volume = (
          getAveragePosition3D(positions, constraint.sites[0]) - delta
        ).dot(
          (
            getAveragePosition3D(positions, constraint.sites[1]) - delta
          ).cross(
            getAveragePosition3D(positions, constraint.sites[2]) - delta
          )
        );

        if( // can this be simplified? -> sign bit XOR?
          ( volume < 0 && constraint.lower > 0)
          || (volume > 0 && constraint.lower < 0)
        ) {
          incorrectNonZeroChiralityConstraints += 1;
        }
      }
    }

    if(nonZeroChiralityConstraints == 0) {
      proportionChiralConstraintsCorrectSign = 1.0;
    }

    proportionChiralConstraintsCorrectSign = static_cast<double>(
      nonZeroChiralityConstraints - incorrectNonZeroChiralityConstraints
    ) / nonZeroChiralityConstraints;

    return proportionChiralConstraintsCorrectSign;
  }

  template<typename Visitor>
  void visitUnfulfilledConstraints(
    const DistanceBoundsMatrix& bounds,
    const Eigen::VectorXd& positions,
    Visitor&& visitor
  ) {
    const unsigned N = positions.size() / dimensionality;

    // Distance bounds deviations
    for(unsigned i = 0; i < N; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        const double ijDistance = (
          positions.template segment<dimensionality>(dimensionality * j)
          - positions.template segment<dimensionality>(dimensionality * i)
        ).norm();

        if(
          ijDistance - bounds.upperBound(i, j) > visitor.deviationThreshold
          || bounds.lowerBound(i, j) - ijDistance > visitor.deviationThreshold
        ) {
          visitor.distanceOverThreshold(i, j, ijDistance);
          if(visitor.earlyExit) {
            return visitor.value;
          }
        }
      }
    }

    // Check chiral bound deviations
    for(const auto& constraint : chiralConstraints) {
      const double volume = (
        getAveragePosition3D(positions, constraint.sites[0])
        - getAveragePosition3D(positions, constraint.sites[3])
      ).dot(
        (
          getAveragePosition3D(positions, constraint.sites[1])
          - getAveragePosition3D(positions, constraint.sites[3])
        ).cross(
          getAveragePosition3D(positions, constraint.sites[2])
          - getAveragePosition3D(positions, constraint.sites[3])
        )
      );

      if(
        volume - constraint.upper > visitor.deviationThreshold
        || constraint.lower - volume > visitor.deviationThreshold
      ) {
        visitor.chiralOverThreshold(constraint, volume);
        if(visitor.earlyExit) {
          return visitor.value;
        }
      }
    }

    // Check dihedral constraints
    for(const auto& constraint : dihedralConstraints) {
      const ThreeDimensionalVector alpha = getAveragePosition3D(positions, constraint.sites[0]);
      const ThreeDimensionalVector beta = getAveragePosition3D(positions, constraint.sites[1]);
      const ThreeDimensionalVector gamma = getAveragePosition3D(positions, constraint.sites[2]);
      const ThreeDimensionalVector delta = getAveragePosition3D(positions, constraint.sites[3]);

      const ThreeDimensionalVector f = alpha - beta;
      const ThreeDimensionalVector g = beta - gamma;
      const ThreeDimensionalVector h = delta - gamma;

      ThreeDimensionalVector a = f.cross(g);
      ThreeDimensionalVector b = h.cross(g);

      double phi = std::atan2(
        a.cross(b).dot(
          -g.normalized()
        ),
        a.dot(b)
      );

      double constraintSumHalved = (constraint.lower + constraint.upper) / 2;

      if(phi < constraintSumHalved - M_PI) {
        phi += 2 * M_PI;
      } else if(phi > constraintSumHalved + M_PI) {
        phi -= 2 * M_PI;
      }

      double term = std::fabs(phi - constraintSumHalved) - (constraint.upper - constraint.lower) / 2;

      if(term > visitor.deviationThreshold) {
        visitor.dihedralOverThreshold(constraint, term);
        if(visitor.earlyExit) {
          return visitor.value;
        }
      }
    }

    return visitor.value;
  }
};

} // namespace DistanceGeometry
} // namespace molassembler
} // namespace Scine

#endif
