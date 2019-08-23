/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Eigen refinement implementation
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_EIGEN_REFINEMENT_PROBLEM_H
#define INCLUDE_MOLASSEMBLER_DG_EIGEN_REFINEMENT_PROBLEM_H

// TMP
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cfenv>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "molassembler/Types.h"
#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"

namespace Scine {
namespace molassembler {
namespace DistanceGeometry {

/**
 * @brief Eigen-based refinement error function
 *
 * @tparam dimensionality 3 or 4 spatial dimensions to refine in
 * @tparam FloatType float or double
 * @tparam SIMD Whether to use rewritten implementations that try to take
 *   advantage of Eigen's SIMD capabilities
 */
template<unsigned dimensionality, typename FloatType, bool SIMD>
class EigenRefinementProblem {
  static_assert(
    dimensionality == 3 || dimensionality == 4,
    "EigenRefinementProblem is only suitable for three or four dimensions"
  );

public:
//!@name Public types
//!@{
  //! Vector layout of positions
  using VectorType = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;

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

//!@name Static member functions
//!@{
  /*! @brief Get a string representation of the type name
   *
   * @complexity{@math{\Theta(1)}}
   */
  static std::string name() {
    std::string name = "Eigen";

    name += "RefinementProblem<dimensionality=";
    name += std::to_string(dimensionality);
    name += ", ";

    if(std::is_same<FloatType, float>::value) {
      name += "float";
    } else {
      name += "double";
    }

    if(SIMD) {
      name += ", SIMD=true>";
    } else {
      name += ", SIMD=false>";
    }

    return name;
  }

  /*! @brief Copy a full dimensional position of an atom
   *
   * @complexity{@math{\Theta(1)}}
   */
  inline static FullDimensionalVector getPosition(
    const VectorType& positions,
    const AtomIndex index
  ) {
    return positions.template segment<dimensionality>(dimensionality * index);
  }

  /*! @brief Copy a three-dimensional position of an atom
   *
   * @complexity{@math{\Theta(1)}}
   */
  inline static ThreeDimensionalVector getPosition3D(
    const VectorType& positions,
    const AtomIndex index
  ) {
    return positions.template segment<3>(dimensionality * index);
  }

  /*! @brief Calculate the average full-dimensional position of several atoms
   *
   * @complexity{@math{\Theta(A)} where A is the number of atoms whose
   * positions to average}
   */
  static FullDimensionalVector getAveragePosition(
    const VectorType& positions,
    const ChiralConstraint::AtomListType& atomList
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

  /*! @brief Calculate the average three-dimensional position of several atoms
   *
   * @complexity{@math{\Theta(A)} where A is the number of atoms whose
   * positions to average}
   */
  static ThreeDimensionalVector getAveragePosition3D(
    const VectorType& positions,
    const ChiralConstraint::AtomListType& atomList
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
  //! Upper distance bounds squared, linearized in i < j
  VectorType upperDistanceBoundsSquared;
  //! Lower distance bounds squared, linearized in i < j
  VectorType lowerDistanceBoundsSquared;
  //! Chiral upper constraints, in sequence of @p chiralConstraints
  VectorType chiralUpperConstraints;
  //! Chiral lower constraints, in sequence of @p chiralConstraints
  VectorType chiralLowerConstraints;
  //! Dihedral bounds' averages, in sequence of @p dihedralConstraints
  VectorType dihedralConstraintSumsHalved;
  //! Dihedral bounds upper minus lower, halved, in sequence
  VectorType dihedralConstraintDiffsHalved;
  //! List of chiral constraints
  std::vector<ChiralConstraint> chiralConstraints;
  //! List of dihedral constraints
  std::vector<DihedralConstraint> dihedralConstraints;
  //! Whether to compress the fourth dimension
  bool compressFourthDimension = false;
  //! Whether to enable dihedral terms
  bool dihedralTerms = false;
//!@}

//!@name Signaling members
//!@{
  mutable double proportionChiralConstraintsCorrectSign = 0.0;
//!@}

//!@name Constructors
//!@{
  EigenRefinementProblem(
    const Eigen::MatrixXd& squaredBounds,
    std::vector<ChiralConstraint> passChiralConstraints,
    std::vector<DihedralConstraint> passDihedralConstraints
  ) : chiralConstraints(std::move(passChiralConstraints)),
      dihedralConstraints(std::move(passDihedralConstraints))
  {
    // Enable floating point exceptions
#ifndef NDEBUG
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    const unsigned N = squaredBounds.cols();
    const unsigned strictlyUpperTriangularElements = N * (N - 1) / 2;

    // Linearize upper distance bounds squared
    upperDistanceBoundsSquared.resize(strictlyUpperTriangularElements);
    lowerDistanceBoundsSquared.resize(strictlyUpperTriangularElements);
    for(unsigned linearIndex = 0, i = 0; i < N; ++i) {
      for(unsigned j = i + 1; j < N; ++j, ++linearIndex) {
        upperDistanceBoundsSquared(linearIndex) = squaredBounds(i, j);
        lowerDistanceBoundsSquared(linearIndex) = squaredBounds(j, i);
      }
    }

    // Vectorize chiral constraint bounds
    const unsigned C = chiralConstraints.size();
    chiralUpperConstraints.resize(C);
    chiralLowerConstraints.resize(C);
    for(unsigned i = 0; i < C; ++i) {
      const ChiralConstraint& constraint = chiralConstraints[i];
      chiralUpperConstraints(i) = constraint.upper;
      chiralLowerConstraints(i) = constraint.lower;
    }

    // Vectorize dihedral constraint bounds sum halves and diff halves
    const unsigned D = dihedralConstraints.size();
    dihedralConstraintSumsHalved.resize(D);
    dihedralConstraintDiffsHalved.resize(D);
    for(unsigned i = 0; i < D; ++i) {
      const DihedralConstraint& constraint = dihedralConstraints[i];
      dihedralConstraintSumsHalved(i) = (constraint.upper + constraint.lower) / 2;
      dihedralConstraintDiffsHalved(i) = (constraint.upper - constraint.lower) / 2;
    }
  }
//!@}

//!@name Contribution functions
//!@{
  /*! @brief Adds pairwise distance error and gradient contributions
   *
   * @complexity{@math{\Omega(N^2)}}
   */
  void distanceContributions(const VectorType& positions, FloatType& error, Eigen::Ref<VectorType> gradient) const {
    // Delegate to SIMD or non-SIMD implementation
    distanceContributionsImpl(positions, error, gradient);
  }

  /*! @brief Adds chiral error and gradient contributions
   *
   * @complexity{@math{\Omega(C)}}
   */
  void chiralContributions(const VectorType& positions, FloatType& error, Eigen::Ref<VectorType> gradient) const {
    // Delegate to SIMD or non-SIMD implementation
    chiralContributionsImpl(positions, error, gradient);
  }

  /*! @brief Adds dihedral error and gradient contributions
   *
   * @complexity{@math{\Omega(D)}}
   */
  void dihedralContributions(const VectorType& positions, FloatType& error, Eigen::Ref<VectorType> gradient) const {
    if(dihedralTerms) {
      // Delegate to SIMD or non-SIMD implementation
      dihedralContributionsImpl(positions, error, gradient);
    }
  }

  /*! @brief Adds fourth dimension error and gradient contributions
   *
   * @complexity{@math{\Theta(N)}}
   */
  void fourthDimensionContributions(const VectorType& positions, FloatType& error, Eigen::Ref<VectorType> gradient) const {
    assert(positions.size() == gradient.size());
    // Collect all fourth dimension values, square and sum, add to error
    if(dimensionality == 4 && compressFourthDimension) {
      const unsigned N = positions.size() / dimensionality;
      for(unsigned i = 0; i < N; ++i) {
        const FloatType fourthDimensionalValue = positions(4 * i + 3);
        error += fourthDimensionalValue * fourthDimensionalValue;
        gradient(4 * i + 3) += 2 * fourthDimensionalValue;
      }
    }
  }
//!@}

  /*!
   * @brief Calculates the error value and gradient for all contributions
   * @param parameters The linearized positions of all particles
   * @pre @p parameters must be evenly divisible by the dimensionality
   * @post Error is stored in @p value and gradient in @p gradient
   */
  void operator() (const VectorType& parameters, FloatType& value, Eigen::Ref<VectorType> gradient) const {
    assert(parameters.size() == gradient.size());
    assert(parameters.size() % dimensionality == 0);

    value = 0;
    gradient.setZero();

    fourthDimensionContributions(parameters, value, gradient);
    dihedralContributions(parameters, value, gradient);
    distanceContributions(parameters, value, gradient);
    chiralContributions(parameters, value, gradient);
  }

  /*! @brief Calculates the number of chiral constraints with correct sign
   *
   * @complexity{@math{\Theta(C)} where @math{C} is the number of chiral
   * constraints}
   */
  double calculateProportionChiralConstraintsCorrectSign(const VectorType& positions) const {
    unsigned nonZeroChiralConstraints = 0;
    unsigned incorrectNonZeroChiralConstraints = 0;

    for(const auto& constraint : chiralConstraints) {
      /* Make sure that chiral constraints meant to cause coplanarity of
       * four indices aren't counted here - Their volume bounds are
       * usually so tight that they might be slightly outside in any given
       * refinement step. If they are outside their target values by a large
       * amount, that is caught by finalStructureAcceptable anyway.
       */
      if(std::fabs(constraint.lower + constraint.upper) > 1e-4) {
        nonZeroChiralConstraints += 1;

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
          incorrectNonZeroChiralConstraints += 1;
        }
      }
    }

    if(nonZeroChiralConstraints == 0) {
      proportionChiralConstraintsCorrectSign = 1.0;
    } else {
      proportionChiralConstraintsCorrectSign = static_cast<double>(
        nonZeroChiralConstraints - incorrectNonZeroChiralConstraints
      ) / nonZeroChiralConstraints;
    }

    return proportionChiralConstraintsCorrectSign;
  }

  /**
   * @brief Visit all unfulfilled constraints
   *
   * @tparam Visitor Type containing a `deviationThreshold` and `result` member
   * variable and implementing the following methods:
   * - distanceOverThreshold: (unsigned, unsigned, double) -> void
   * - chiralOverThreshold: (ChiralConstraint, double) -> void
   * - dihedralOverThreshold: (DihedralConstraint, double) -> void
   * @param bounds distance bounds instance
   * @param positions parameter positions to the error function
   * @param visitor Instance of VisitorType
   *
   * This function is very helpful in implementing multiple behaviors, like
   * debug or error information on unsatisfied bonds or checking whether all
   * constraints are within their bounds.
   *
   * @returns `visitor.result`
   */
  template<typename Visitor>
  auto visitUnfulfilledConstraints(
    const DistanceBoundsMatrix& bounds,
    const VectorType& positions,
    Visitor&& visitor
  ) const {
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

      const ThreeDimensionalVector a = f.cross(g);
      const ThreeDimensionalVector b = h.cross(g);

      const double dihedral = std::atan2(
        a.cross(b).dot(
          -g.normalized()
        ),
        a.dot(b)
      );

      const double constraintSumHalved = (constraint.lower + constraint.upper) / 2;

      double phi;
      if(dihedral < constraintSumHalved - M_PI) {
        phi = dihedral + 2 * M_PI;
      } else if(dihedral > constraintSumHalved + M_PI) {
        phi = dihedral - 2 * M_PI;
      } else {
        phi = dihedral;
      }

      const double term = std::fabs(phi - constraintSumHalved) - (constraint.upper - constraint.lower) / 2;

      if(term > visitor.deviationThreshold) {
        visitor.dihedralOverThreshold(constraint, term);
        if(visitor.earlyExit) {
          return visitor.value;
        }
      }
    }

    return visitor.value;
  }

private:
//!@name Contribution implementations
//!@{
  /*!
   * @brief Adds distance error and gradient contributions (non-SIMD)
   */
  template<bool dependent = SIMD, std::enable_if_t<!dependent, int>...>
  void distanceContributionsImpl(
    const VectorType& positions,
    FloatType& error,
    Eigen::Ref<VectorType> gradient
  ) const {
    assert(positions.size() == gradient.size());
    const unsigned N = positions.size() / dimensionality;

    for(unsigned linearIndex = 0, i = 0; i < N - 1; ++i) {
      for(unsigned j = i + 1; j < N; ++j, ++linearIndex) {
        const FloatType lowerBoundSquared = lowerDistanceBoundsSquared(linearIndex);
        const FloatType upperBoundSquared = upperDistanceBoundsSquared(linearIndex);
        assert(lowerBoundSquared <= upperBoundSquared);

        // For both
        const FullDimensionalVector positionDifference = (
          positions.template segment<dimensionality>(dimensionality * i)
          - positions.template segment<dimensionality>(dimensionality * j)
        );

        const FloatType squareDistance = positionDifference.squaredNorm();

        // Upper term
        const FloatType upperTerm = squareDistance / upperBoundSquared - 1;

        if(upperTerm > 0) {
          error += upperTerm * upperTerm;

          const FullDimensionalVector f = 4 * positionDifference * upperTerm / upperBoundSquared;

          gradient.template segment<dimensionality>(dimensionality * i) += f;
          gradient.template segment<dimensionality>(dimensionality * j) -= f;
        } else {
          // Lower term is only possible if the upper term does not contribute
          const FloatType quotient = lowerBoundSquared + squareDistance;
          const FloatType lowerTerm = 2 * lowerBoundSquared / quotient - 1;

          if(lowerTerm > 0) {
            error += lowerTerm * lowerTerm;

            const FullDimensionalVector g = 8 * lowerBoundSquared * positionDifference * lowerTerm / (
              quotient * quotient
            );

            /* We use -= because the lower term needs the position vector
             * difference (j - i), so we reuse positionDifference and just subtract
             * from the gradient instead of adding to it
             */
            gradient.template segment<dimensionality>(dimensionality * i) -= g;
            gradient.template segment<dimensionality>(dimensionality * j) += g;
          }
        }
      }
    }
  }

  /*!
   * @brief SIMD implementation of distance contributions
   */
  template<bool dependent = SIMD, std::enable_if_t<dependent, int>...>
  void distanceContributionsImpl(
    const VectorType& positions,
    FloatType& error,
    Eigen::Ref<VectorType> gradient
  ) const {
    assert(positions.size() == gradient.size());
    /* Using the squared distance bounds to avoid unneeded square-root calculations
     *
     * NOTE: positions are full-dimensional
     *
     * for each nonredundant pair of atoms (i < j)
     *   Calculate a squared vector distance:
     *   sq_distance = (position[j] - position[i]).array().square().sum()
     *
     *   upper_term = (sq_distance / upperBoundSquared) - 1
     *
     *   if upper_term > 0:
     *     error += upper_term.square()
     *
     *   lower_term = (
     *    2 * lowerBoundSquared / (lowerBoundSquared + sq_distance)
     *   ) - 1
     *
     *   if lowerTerm > 0:
     *     error += lower_term.square()
     *
     * SIMDize:
     *   sq_distances = list of upper triangular square distance terms
     *   upper_bounds_sq = list of corresponding upper bounds (identical across runs, no need to gather multiple times)
     *   lower_bounds_sq = list of corresponding lower bounds (same as above)
     *
     *   upper_terms = (sq_distances / upper_bounds_sq) - 1
     *   error += upper_terms.max(0.0).square().sum()
     *
     *   lower_terms = (2 * lower_bounds_sq / (lower_bounds_sq + sq_distances)) - 1
     *   error += lower_terms.max(0.0).square().sum()
     *
     */

    // Calculate all full-dimensional N(N-1) / 2 position differences
    const unsigned N = positions.size() / dimensionality;
    const unsigned differencesCount = N * (N - 1) / 2;

    FullDimensionalMatrixType positionDifferences(dimensionality, differencesCount);
    {
      unsigned offset = 0;
      for(unsigned i = 0; i < N - 1; ++i) {
        auto iPosition = positions.template segment<dimensionality>(dimensionality * i);
        for(unsigned j = i + 1; j < N; ++j) {
          positionDifferences.col(offset) = (
            iPosition
            - positions.template segment<dimensionality>(dimensionality * j)
          );
          ++offset;
        }
      }
    }

    // SIMD
    const VectorType squareDistances = positionDifferences.colwise().squaredNorm();

    // SIMD
    const VectorType upperTerms = squareDistances.array() / upperDistanceBoundsSquared.array() - 1;

#ifndef NDEBUG
    {
      // Check correctness of results so far
      unsigned offset = 0;
      for(unsigned i = 0; i < N - 1; ++i) {
        for(unsigned j = i + 1; j < N; ++j) {
          const FullDimensionalVector positionDifference = (
            positions.template segment<dimensionality>(dimensionality * i)
            - positions.template segment<dimensionality>(dimensionality * j)
          );

          assert(positionDifferences.col(offset) == positionDifference);

          const FloatType squareDistance = positionDifference.squaredNorm();

          assert(squareDistance == squareDistances(offset));

          // Upper term
          const FloatType upperTerm = squareDistance / upperDistanceBoundsSquared(offset) - 1;

          assert(upperTerm == upperTerms(offset));

          ++offset;
        }
      }
    }
#endif

    // Traverse upper terms and calculate gradients
    for(unsigned iOffset = 0, i = 0; i < N - 1; ++i) {
      const unsigned crossTerms = N - i - 1;
      for(unsigned jOffset = 0, j = i + 1;  j < N; ++j, ++jOffset) {
        const FloatType upperTerm = upperTerms(iOffset + jOffset);
        if(upperTerm > 0) {
          error += upperTerm * upperTerm;

          /* This i-j combination has an upper value contribution
           *
           * calculate gradient and apply to i and j
           */
          const FullDimensionalVector f = (
            4 * positionDifferences.col(iOffset + jOffset) * upperTerm
            / upperDistanceBoundsSquared(iOffset + jOffset)
          );

          // position difference is i - j, not j - i (!)
          gradient.template segment<dimensionality>(dimensionality * i) += f;
          gradient.template segment<dimensionality>(dimensionality * j) -= f;
        } else {
          /* This i-j combination MAYBE has a lower value contribution
           */
          const auto& lowerBoundSquared = lowerDistanceBoundsSquared(iOffset + jOffset);
          const FloatType quotient = lowerBoundSquared + squareDistances(iOffset + jOffset);
          const FloatType lowerTerm = 2 * lowerBoundSquared / quotient - 1;

          if(lowerTerm > 0) {
            // Add it to the error
            error += lowerTerm * lowerTerm;

            // Calculate gradient and apply
            const FullDimensionalVector g = (
              8 * lowerBoundSquared * positionDifferences.col(iOffset + jOffset) * lowerTerm / (
                quotient * quotient
              )
            );

            // position difference is i - j, not j - i (!)
            gradient.template segment<dimensionality>(dimensionality * i) -= g;
            gradient.template segment<dimensionality>(dimensionality * j) += g;
          }
        }
      }

      iOffset += crossTerms;
    }
  }

  /*!
   * @brief Adds chiral error and gradient contributions
   */
  template<bool dependent = SIMD, std::enable_if_t<!dependent, int>...>
  void chiralContributionsImpl(const VectorType& positions, FloatType& error, Eigen::Ref<VectorType> gradient) const {
    unsigned nonZeroChiralConstraints = 0;
    unsigned incorrectNonZeroChiralConstraints = 0;

    for(const auto& constraint : chiralConstraints) {
      const ThreeDimensionalVector alpha = getAveragePosition3D(positions, constraint.sites[0]);
      const ThreeDimensionalVector beta = getAveragePosition3D(positions, constraint.sites[1]);
      const ThreeDimensionalVector gamma = getAveragePosition3D(positions, constraint.sites[2]);
      const ThreeDimensionalVector delta = getAveragePosition3D(positions, constraint.sites[3]);

      const ThreeDimensionalVector alphaMinusDelta = alpha - delta;
      const ThreeDimensionalVector betaMinusDelta = beta - delta;
      const ThreeDimensionalVector gammaMinusDelta = gamma - delta;

      const FloatType volume = (alphaMinusDelta).dot(
        betaMinusDelta.cross(gammaMinusDelta)
      );

      if(std::fabs(constraint.lower + constraint.upper) > 1e-4) {
        ++nonZeroChiralConstraints;
        if(
          ( volume < 0 && constraint.lower > 0)
          || (volume > 0 && constraint.lower < 0)
        ) {
          ++incorrectNonZeroChiralConstraints;
        }
      }

      // Upper bound term
      const FloatType upperTerm = volume - static_cast<FloatType>(constraint.upper);

      if(upperTerm > 0) {
        error += upperTerm * upperTerm;
      }

      // Lower bound term
      const FloatType lowerTerm = static_cast<FloatType>(constraint.lower) - volume;
      if(lowerTerm > 0) {
        error += lowerTerm * lowerTerm;
      }

      const FloatType factor = 2 * (
        std::max(FloatType {0}, upperTerm)
        - std::max(FloatType {0}, lowerTerm)
      );

      /* Make sure computing all cross products is worth it by checking that
       * one of both terms is actually greater than 0
       */
      if(factor != 0) {
        // Specific to deltaI only but still repeated
        const ThreeDimensionalVector alphaMinusGamma = alpha - gamma;

        const ThreeDimensionalVector betaMinusGamma = beta - gamma;

        const ThreeDimensionalVector iContribution = (
          factor * betaMinusDelta.cross(gammaMinusDelta)
          / constraint.sites[0].size()
        );

        const ThreeDimensionalVector jContribution = (
          factor * gammaMinusDelta.cross(alphaMinusDelta)
          / constraint.sites[1].size()
        );

        const ThreeDimensionalVector kContribution = (
          factor * alphaMinusDelta.cross(betaMinusDelta)
          / constraint.sites[2].size()
        );

        const ThreeDimensionalVector lContribution = (
          factor * betaMinusGamma.cross(alphaMinusGamma)
          / constraint.sites[3].size()
        );

        for(const AtomIndex alphaI : constraint.sites[0]) {
          gradient.template segment<3>(dimensionality * alphaI) += iContribution;
        }

        for(const AtomIndex betaI : constraint.sites[1]) {
          gradient.template segment<3>(dimensionality * betaI) += jContribution;
        }

        for(const AtomIndex gammaI : constraint.sites[2]) {
          gradient.template segment<3>(dimensionality * gammaI) += kContribution;
        }

        for(const AtomIndex deltaI : constraint.sites[3]) {
          gradient.template segment<3>(dimensionality * deltaI) += lContribution;
        }
      }
    }

    // Set signaling member
    if(nonZeroChiralConstraints == 0) {
      proportionChiralConstraintsCorrectSign = 1;
    } else {
      proportionChiralConstraintsCorrectSign = static_cast<double>(
        nonZeroChiralConstraints - incorrectNonZeroChiralConstraints
      ) / nonZeroChiralConstraints;
    }
  }

  /*!
   * @brief SIMD implementation of chiral contributions
   */
  template<bool dependent = SIMD, std::enable_if_t<dependent, int>...>
  void chiralContributionsImpl(const VectorType& positions, FloatType& error, Eigen::Ref<VectorType> gradient) const {
    assert(positions.size() == gradient.size());
    /* We have to calculate the adjusted signed volume for each chiral
     * constraint once for both value and gradient, so we could principally do
     * all of them at once
     *
     * Each site sequence has four averaged three-dimensional positions as
     * source data, need to do three vector differences:
     *
     * minuend = sites[0], sites[1], sites[2]
     * subtrahend = sites[3], sites[3], sites[3]
     * difference = minuend - subtrahend (SIMD advantage the more of these we do at once)
     *
     * volume = difference[0].dot(difference[1].cross(difference[2]));
     *
     * then, we do:
     *   upperTerm = volume - constraint.upper
     *   if(upperTerm > 0) {
     *     error += upperTerm * upperTerm;
     *   }
     *
     * in SIMD baby steps:
     *   volumes = volume of constraint 0, ...
     *   upper_constraints = constraint.upper of constraint 0, ...
     *   upper_term = volumes - upper_constraints
     *   positive_upper_terms = upper_term.cwiseMax(0.0)
     *   positive_upper_terms_sq_arr = positive_upper_terms.array().square()
     *   positive_upper_term_sq_sum = positive_upper_terms_sq_arr.sum()
     *   error += positive_upper_term_sq_sum
     *
     *   or altogether:
     *     error += (volumes - upper_constraints).array().max(0.0).square().sum()
     *
     * and similarly for the lower term:
     *   error += (lower_constraints - volumes).array().max(0.0).square().sum()
     *
     * gradient:
     * calculate adjusted signed volume of constraint (see chiralError)
     * calculate upper and lower term (see chiralError)
     *
     * if either upper or lower term are greater than zero,
     *   factor = 2 * (
     *     max(0.0, upperTerm) - max(0.0, lowerTerm)
     *   )
     *
     * SIMD-izing this if is a little difficult. but calculating factor for all
     * chiral constraints would be
     *   factor = (upper_terms.max(0.0) - lower_terms.max(0.0)) * 2
     *
     * Then collect all chiral constraints for which factor != 0.0 and for each
     * of those, we calculate:
     *
     *   a) alpha - delta = sites[0] - sites[3]
     *   b) beta - delta = sites[1] - sites[3]
     *   c) gamma - delta = sites[2] - sites[3]
     *   d) alpha - gamma = sites[0] - sites[2]
     *   e) beta - gamma = sites[1] - sites[2]
     *
     *   can SIMD-ize this part with
     *     minuend = sites[0], sites[1], sites[2], sites[0], sites[1] (gather)
     *     subtrahend = sites[3], sites[3], sites[3], sites[2], sites[2] (gather)
     *     difference = minuend - subtrahend (SIMD op)
     *     a = difference[0], b = difference[1], ...  (scatter)
     *
     *   ijkl contributions
     *     individually: factor * (cross product of two a-e terms) / site size count (different for each term)
     *
     *     cross_products (arr) = b.cross(c), c.cross(a), a.cross(b), e.cross(d)
     *     site_sizes (arr) = sites[0].size(), sites[1].size(), sites[2].size(), sites[3].size()
     *     ijkl_contributions = (crossProducts / site_sizes) * factor
     *
     *   then distribute each contribution vector across all component atoms of each
     *   site (i to site 0, etc.)
     *
     * this can (in principle) be done for all chiral constraints at once,
     * except for the distribution step, which then must be done serially
     */
    const unsigned C = chiralConstraints.size();
    VectorType volumes(C);

    ThreeDimensionalMatrixType positionDifferences, minuend;
    {
      // Three differences per constraint * number of constraints
      unsigned differencesCount = 3 * chiralConstraints.size();
      ThreeDimensionalMatrixType subtrahend;
      minuend.resize(Eigen::NoChange, differencesCount);
      subtrahend.resize(Eigen::NoChange, differencesCount);

      // Populate minuend and subtrahend
      for(unsigned constraintIndex = 0; constraintIndex < C; ++constraintIndex) {
        const ChiralConstraint& constraint = chiralConstraints[constraintIndex];

        auto lastSitePosition = getAveragePosition3D(positions, constraint.sites.back());
        for(unsigned siteIndex = 0; siteIndex < 3; ++siteIndex) {
          const unsigned offset = constraintIndex * 3 + siteIndex;
          minuend.col(offset) = getAveragePosition3D(
            positions,
            constraint.sites[siteIndex]
          );
          subtrahend.col(offset) = lastSitePosition;
        }
      }

      // SIMD (benchmark, of questionable value!)
      positionDifferences = minuend - subtrahend;
    } // release subtrahend

    // Reduce positionDifferences to volumes
    for(unsigned constraintIndex = 0; constraintIndex < C; ++constraintIndex) {
      volumes(constraintIndex) = positionDifferences.col(
        constraintIndex * 3 + 0
      ).dot(
        positionDifferences.col(constraintIndex * 3 + 1).cross(
          positionDifferences.col(constraintIndex * 3 + 2)
        )
      );
    }

    // Count constraints with correct signs and set signaling member
    unsigned nonZeroChiralConstraints = 0;
    unsigned incorrectNonZeroChiralConstraints = 0;
    for(unsigned constraintIndex = 0; constraintIndex < C; ++constraintIndex) {
      const ChiralConstraint& constraint = chiralConstraints[constraintIndex];
      if(std::fabs(constraint.lower + constraint.upper) > 1e-4) {
        ++nonZeroChiralConstraints;
        const double volume = volumes(constraintIndex);
        if(
          ( volume < 0 && constraint.lower > 0)
          || (volume > 0 && constraint.lower < 0)
        ) {
          ++incorrectNonZeroChiralConstraints;
        }
      }
    }
    if(nonZeroChiralConstraints == 0) {
      proportionChiralConstraintsCorrectSign = 1;
    } else {
      proportionChiralConstraintsCorrectSign = static_cast<double>(
        nonZeroChiralConstraints - incorrectNonZeroChiralConstraints
      ) / nonZeroChiralConstraints;
    }

    // SIMD
    const VectorType upperTerms = volumes - chiralUpperConstraints;
    error += upperTerms.array().max(0.0).square().sum();

    // SIMD
    const VectorType lowerTerms = chiralLowerConstraints - volumes;
    error += lowerTerms.array().max(0.0).square().sum();

    // SIMD
    const VectorType factors = 2 * (
      upperTerms.array().max(0.0) - lowerTerms.array().max(0.0)
    ).matrix();

    /* A priori unknown how many chiral constraints contribute to the error
     * function and gradients over time and how much value an optimization would
     * have:
     *
     * Possible optimization: collect all constraints that have factor != 0.0
     * then populate ijklContributions and site_sizes
     */
    ThreeDimensionalMatrixType ijklContributions(3, 4 * C);
    VectorType siteSizes(4 * C);
    for(unsigned constraintIndex = 0; constraintIndex < C; ++constraintIndex) {
      const FloatType factor = factors(constraintIndex);
      if(factor != 0) {
        // Collect sub-expressions for the various terms, does nothing yet
        auto a = positionDifferences.col(3 * constraintIndex + 0);
        auto b = positionDifferences.col(3 * constraintIndex + 1);
        auto c = positionDifferences.col(3 * constraintIndex + 2);
        auto d = (
          minuend.col(3 * constraintIndex + 0)
          - minuend.col(3 * constraintIndex + 2)
        );
        auto e = (
          minuend.col(3 * constraintIndex + 1)
          - minuend.col(3 * constraintIndex + 2)
        );

        // ijkl contributions
        ijklContributions.col(4 * constraintIndex + 0) = b.cross(c);
        ijklContributions.col(4 * constraintIndex + 1) = c.cross(a);
        ijklContributions.col(4 * constraintIndex + 2) = a.cross(b);
        ijklContributions.col(4 * constraintIndex + 3) = e.cross(d);

        const ChiralConstraint& constraint = chiralConstraints[constraintIndex];

        // site sizes
        siteSizes(4 * constraintIndex + 0) = constraint.sites[0].size();
        siteSizes(4 * constraintIndex + 1) = constraint.sites[1].size();
        siteSizes(4 * constraintIndex + 2) = constraint.sites[2].size();
        siteSizes(4 * constraintIndex + 3) = constraint.sites[3].size();
      }
    }

    // SIMD
    /* - ijklContributions is 3 x (4 * C)
     * - factors is 1xC
     * - siteSizes is 1 x (4 * C)
     *
     * -> for each constraint, multiply the corresponding 3x4 ijkl block
     *    by (factors * siteSizes).asDiagonal()
     */
    for(unsigned constraintIndex = 0; constraintIndex < C; ++constraintIndex) {
      // SIMD (?)
      ijklContributions.template block<3, 4>(0, 4 * constraintIndex) *= (
        siteSizes.template segment<4>(4 * constraintIndex) * factors(constraintIndex)
      ).asDiagonal();
    }
    // ijklContributions = ijklContributions * (factors.array() * siteSizes.array()).matrix().asDiagonal();
    // ijklContributions.array().colwise() *= (factors.array() * siteSizes.array()).transpose();

    // Distribute ijkl contributions into gradient
    for(unsigned constraintIndex = 0; constraintIndex < C; ++constraintIndex) {
      if(factors(constraintIndex) != 0) {
        const ChiralConstraint& constraint = chiralConstraints[constraintIndex];
        /* ijklContributions.col is a Eigen::ColExpr.
         * If FloatType == double, then iContribution is Eigen::ColExpr -> no copy
         * If FloatType == float, then iContribution is Eigen::ColExpr.cast<double>().eval() == Vector3d -> copying cast
         */
        auto iContribution = ijklContributions.col(4 * constraintIndex + 0);
        auto jContribution = ijklContributions.col(4 * constraintIndex + 1);
        auto kContribution = ijklContributions.col(4 * constraintIndex + 2);
        auto lContribution = ijklContributions.col(4 * constraintIndex + 3);

        for(const AtomIndex alphaI : constraint.sites[0]) {
          gradient.template segment<3>(dimensionality * alphaI) += iContribution;
        }

        for(const AtomIndex betaI : constraint.sites[1]) {
          gradient.template segment<3>(dimensionality * betaI) += jContribution;
        }

        for(const AtomIndex gammaI : constraint.sites[2]) {
          gradient.template segment<3>(dimensionality * gammaI) += kContribution;
        }

        for(const AtomIndex deltaI : constraint.sites[3]) {
          gradient.template segment<3>(dimensionality * deltaI) += lContribution;
        }
      }
    }
  }

  /*!
   * @brief Adds dihedral error and gradient contributions
   */
  template<bool dependent = SIMD, std::enable_if_t<!dependent, int>...>
  void dihedralContributionsImpl(const VectorType& positions, FloatType& error, Eigen::Ref<VectorType> gradient) const {
    assert(positions.size() == gradient.size());
    constexpr FloatType reductionFactor = 1.0 / 10;

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

      // Calculate the dihedral angle
      FloatType phi = std::atan2(
        a.cross(b).dot(
          -g.normalized()
        ),
        a.dot(b)
      );

      const FloatType constraintSumHalved = (static_cast<FloatType>(constraint.lower) + static_cast<FloatType>(constraint.upper)) / 2;

      constexpr FloatType pi {M_PI};

      if(phi < constraintSumHalved - pi) {
        phi += 2 * pi;
      } else if(phi > constraintSumHalved + pi) {
        phi -= 2 * pi;
      }

      const FloatType w_phi = phi - constraintSumHalved;

      FloatType h_phi = std::fabs(w_phi) - (static_cast<FloatType>(constraint.upper) - static_cast<FloatType>(constraint.lower)) / 2;

      /* "Apply the max function": If h <= 0, then the max function yields zero,
       * so we can skip this constraint
       */
      if(h_phi <= 0) {
        continue;
      }

      // Error contribution
      error += h_phi * h_phi * reductionFactor;

      // Multiply with sgn (w)
      h_phi *= static_cast<int>(0 < w_phi) - static_cast<int>(w_phi < 0);

      // Multiply in 2 and the reduction factor
      h_phi *= 2 * reductionFactor;

      // Precompute some reused expressions
      const FloatType gLength = g.norm();
      if(gLength == 0) {
        throw std::runtime_error("Encountered zero-length g vector in dihedral errors");
      }

      const FloatType aLengthSq = a.squaredNorm();
      if(aLengthSq == 0) {
        throw std::runtime_error("Encountered zero-length a vector in dihedral gradient contributions");
      }
      a /= aLengthSq;
      const FloatType bLengthSq = b.squaredNorm();
      if(bLengthSq == 0) {
        throw std::runtime_error("Encountered zero-length b vector in dihedral gradient contributions");
      }
      b /= bLengthSq;
      const FloatType fDotG = f.dot(g);
      const FloatType gDotH = g.dot(h);

      const ThreeDimensionalVector iContribution = (
        -(h_phi / constraint.sites[0].size()) * gLength * a
      );

      const ThreeDimensionalVector jContribution = (
        (h_phi / constraint.sites[1].size()) * (
          (gLength + fDotG / gLength) * a
          - (gDotH / gLength) * b
        )
      );

      const ThreeDimensionalVector kContribution = (
        (h_phi / constraint.sites[2].size()) * (
          (gDotH / gLength - gLength) * b
          - (fDotG / gLength) * a
        )
      );

      const ThreeDimensionalVector lContribution = (
        (h_phi / constraint.sites[3].size()) * gLength * b
      );

      if(
        !iContribution.array().isFinite().all()
        || !jContribution.array().isFinite().all()
        || !kContribution.array().isFinite().all()
        || !lContribution.array().isFinite().all()
      ) {
        throw std::runtime_error("Encountered non-finite dihedral gradient contributions");
      }

      /* Distribute contributions among constituting indices */
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
   * @brief SIMD implementation of dihedral contributions
   */
  template<bool dependent = SIMD, std::enable_if_t<dependent, int>...>
  void dihedralContributionsImpl(const VectorType& positions, FloatType& error, Eigen::Ref<VectorType> gradient) const {
    assert(positions.size() == gradient.size());
    /* NOTE: site is average position of constituting atoms in three dimensions
     *
     * For each dihedral constraint,
     *   calculate the current dihedral angle (very little SIMD value):
     *   f = site[1] - site[0]
     *   g = site[2] - site[1]
     *   h = site[3] - site[2]
     *   a = f.cross(g)
     *   b = g.cross(h)
     *   phi = atan2(
     *     a.cross(b).dot(
     *       g.normalized()
     *     ),
     *     a.dot(b)
     *   )
     *
     *   constraint_sum_halved = (constraint.upper + constraint.lower) / 2
     *   constraint_diff_halved = (constraint.upper - constraint.lower) / 2
     *
     *   if phi < constraint_sum_halved - pi:
     *     phi += 2 * pi
     *   elif phi > constraint_sum_halved + pi:
     *     phi -= 2 * pi
     *
     *   term = fabs(phi - constraint_sum_halved) - constraint_diff_halved
     *
     *   if term > 0:
     *     error += term * term / 10
     *
     * SIMD:
     *   use constraint_sums_halved and constraint_diffs_halved, calculated at ctor
     *
     *   phis = calculate dihedral and adjust by 2 pi (perhaps with select?)
     *
     *   adjust phi by 2pi as before
     *
     *   w_phi = phi - constraint_sum_halved
     *   h_phi = fabs(w_phi) - constraint_diff_halved
     *   error += h_phi.max(0.0).square().sum() / 10 (scale reduction factor)
     *
     * NOTE
     *   When reading this, remember that matrix.col(i) yields a ColExpr, a
     *   proxy object instance that directly references matrix' data without a
     *   copy. The use of these objects manipulate the matrix' underlying data!
     */
    constexpr FloatType reductionFactor = 1.0 / 10;
    const unsigned D = dihedralConstraints.size();
    ThreeDimensionalMatrixType fs(3, D), gs(3, D), hs(3, D), as(3, D), bs(3, D);
    VectorType phis(D);
    for(unsigned i = 0; i < D; ++i) {
      const DihedralConstraint& constraint = dihedralConstraints[i];

      const ThreeDimensionalVector alpha = getAveragePosition3D(positions, constraint.sites[0]);
      const ThreeDimensionalVector beta = getAveragePosition3D(positions, constraint.sites[1]);
      const ThreeDimensionalVector gamma = getAveragePosition3D(positions, constraint.sites[2]);
      const ThreeDimensionalVector delta = getAveragePosition3D(positions, constraint.sites[3]);

      // Get proxy objects for column vectors in the matrices
      auto f = fs.col(i);
      auto g = gs.col(i);
      auto h = hs.col(i);
      auto a = as.col(i);
      auto b = bs.col(i);

      /* Each of these expressions assigns a column in one of the two-character
       * matrices
       */
      f = alpha - beta;
      g = beta - gamma;
      h = delta - gamma;

      // WARNING: a and b can be zero-length-vectors!
      a = f.cross(g);
      b = h.cross(g);

      // Calculate the dihedral angle
      FloatType& phi = phis(i);
      phi = std::atan2(
        a.cross(b).dot(
          -g.normalized()
        ),
        a.dot(b)
      );

      constexpr FloatType pi {M_PI};

      if(phi < dihedralConstraintSumsHalved(i) - pi) {
        phi += 2 * pi;
      } else if(phi > dihedralConstraintSumsHalved(i) + pi) {
        phi -= 2 * pi;
      }
    }

    // SIMD
    const VectorType w_phis = phis - dihedralConstraintSumsHalved;

    // SIMD
    VectorType h_phis = w_phis.cwiseAbs() - dihedralConstraintDiffsHalved;

    // SIMD
    error += h_phis.array().max(0.0).square().sum() / 10;

    // Gradient contribution
    for(unsigned i = 0; i < D; ++i) {
      FloatType& h_phi = h_phis(i);

      // Early exit conditions (max function yields zero)
      if(h_phi <= 0) {
        continue;
      }

      /* Skip this constraint if length of g is zero -> all contributions would
       * be NaNs since we divide by it.
       */
      const ThreeDimensionalVector& g = gs.col(i);
      const FloatType gLength = g.norm();
      if(gLength == 0) {
        throw std::runtime_error("Encountered zero-length g vector in dihedral errors");
      }

      const FloatType& w_phi = w_phis(i);
      // Multiply with sgn (w)
      h_phi *= static_cast<int>(0 < w_phi) - static_cast<int>(w_phi < 0);

      // Multiply in 2 and the reduction factor
      h_phi *= 2 * reductionFactor;

      const DihedralConstraint& constraint = dihedralConstraints[i];
      const ThreeDimensionalVector& f = fs.col(i);
      const ThreeDimensionalVector& h = hs.col(i);
      auto a = as.col(i);
      auto b = bs.col(i);

      const FloatType aLengthSquared = a.squaredNorm();
      if(aLengthSquared == 0) {
        throw std::runtime_error("Encountered zero-length a vector in dihedral gradient contributions");
      }
      a /= aLengthSquared;

      const FloatType bLengthSquared = b.squaredNorm();
      if(bLengthSquared == 0) {
        throw std::runtime_error("Encountered zero-length b vector in dihedral gradient contributions");
      }
      b /= bLengthSquared;

      const FloatType fDotG = f.dot(g);
      const FloatType gDotH = g.dot(h);

      const ThreeDimensionalVector iContribution = (
        -(h_phi / constraint.sites[0].size()) * gLength * a
      );

      const ThreeDimensionalVector jContribution = (
        (h_phi / constraint.sites[1].size()) * (
          (gLength + fDotG / gLength) * a
          - (gDotH / gLength) * b
        )
      );

      const ThreeDimensionalVector kContribution = (
        (h_phi / constraint.sites[2].size()) * (
          (gDotH / gLength - gLength) * b
          - (fDotG / gLength) * a
        )
      );

      const ThreeDimensionalVector lContribution = (
        (h_phi / constraint.sites[3].size()) * gLength * b
      );

      if(
        !iContribution.array().isFinite().all()
        || !jContribution.array().isFinite().all()
        || !kContribution.array().isFinite().all()
        || !lContribution.array().isFinite().all()
      ) {
        throw std::runtime_error("Encountered non-finite dihedral gradient contributions");
      }

      /* Distribute contributions among constituting indices */
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
//!@}
};

/**
 * @brief This is for when you have a fully qualified Refinement typename but
 *   want to get the template arguments back
 */
template<typename RefinementType>
struct RefinementTraits {
  // Trick to make the compiler tell me the template arguments to RefinementType
  template<
    template<unsigned, typename, bool> class BaseClass,
    unsigned dimensionality,
    typename FloatingPointType,
    bool SIMD
  > static auto templateArgumentHelper(BaseClass<dimensionality, FloatingPointType, SIMD> /* a */)
    -> std::tuple<
      std::integral_constant<unsigned, dimensionality>,
      FloatingPointType,
      std::integral_constant<bool, SIMD>
    >
  {
    return {};
  }

  using TemplateArgumentsTuple = decltype(templateArgumentHelper(std::declval<RefinementType>()));
  using DimensionalityConstant = std::tuple_element_t<0, TemplateArgumentsTuple>;
  using FloatingPointType = std::tuple_element_t<1, TemplateArgumentsTuple>;
  using SIMDConstant = std::tuple_element_t<2, TemplateArgumentsTuple>;
};

} // namespace DistanceGeometry
} // namespace molassembler
} // namespace Scine

#endif
