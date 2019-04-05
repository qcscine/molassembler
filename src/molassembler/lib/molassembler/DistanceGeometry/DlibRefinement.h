/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Specifies the refinement minimization problem for the dlib library.
 */

#ifndef INCLUDE_MOLASSEMBLER_DG_DLIB_REFINEMENT_PROBLEM_H
#define INCLUDE_MOLASSEMBLER_DG_DLIB_REFINEMENT_PROBLEM_H

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"

#include <dlib/optimization.h>

namespace Scine {
namespace molassembler {
namespace DistanceGeometry {

/**
 * @brief Functor to calculate the DG refinement error function value
 *
 * Refinement progresses in two stages:
 * - Uncompressed, in which the spatial coordinates are free to expand into
 *   the fourth dimension. This is so that chiral constraints can invert
 *   more easily. Once all chiral constraints have the correct sign, we move
 *   on to the next stage.
 * - Compressed, in which the remaining distance bounds and chiral constraints
 *   are optimized further and the fourth dimension is "compressed" out.
 */
struct ErrorFunctionValue {
  // Typedefs
  using Vector = dlib::matrix<double, 0, 1>;

//!@name Static member functions
//!@{
  /**
   * @brief Fetches a three-dimensional position sub-part of the positions vector
   *
   * @param positions All positions
   * @param i The index of the particle whose position to extract
   *
   * @return The position of a particle in the positions
   */
  inline static dlib::vector<double, 3> getPos3D(const Vector& positions, const unsigned i) {
    assert(4 * i + 3 < positions.size());

    return dlib::rowm(
      positions,
      dlib::range(4 * i, 4 * i + 2)
    );
  }

  /**
   * @brief Fetches a four-dimensional position sub-part of the positions vector
   *
   * @param positions The long vector of positions from which to extract a
   *   four-dimensional position
   * @param i The index of the particle whose position to extract
   *
   * @return The four-dimensional position of a particle in the positions
   */
  inline static Vector getPos(const Vector& positions, const unsigned i) {
    assert(4 * i + 3 < positions.size());

    return dlib::rowm(
      positions,
      dlib::range(4 * i, 4 * i + 3)
    );
  }

  /**
   * @brief Calculates the average three-dimensional position of a set of particles
   *
   * @param positions All positions
   * @param atomList A list of atom indices whose average position to calculate
   *
   * @return The average three-dimensional position of a set of particles
   */
  inline static dlib::vector<double, 3> getAveragePos3D(
    const Vector& positions,
    const ChiralityConstraint::AtomListType& atomList
  ) {
    if(atomList.size() == 1) {
      return getPos3D(positions, atomList.front());
    }

    dlib::vector<double, 3> sum {0.0, 0.0, 0.0};

    for(const auto& atomIndex : atomList) {
      sum += dlib::rowm(positions, dlib::range(4 * atomIndex, 4 * atomIndex + 2));
    }

    sum /= atomList.size();

    return sum;
  }

  /**
   * @brief Calculates the average four-dimensional position of a set of particles
   *
   * @param positions All positions
   * @param atomList A list of atom indices whose average position to calculate
   *
   * @return The average four-dimensional position of a set of particles
   */
  inline static Vector getAveragePos(
    const Vector& positions,
    const ChiralityConstraint::AtomListType& atomList
  ) {
    if(atomList.size() == 1) {
      return getPos(positions, atomList.front());
    }

    Vector sum = dlib::zeros_matrix<double>(4, 1);

    for(const auto& atomIndex : atomList) {
      sum += dlib::rowm(positions, dlib::range(4 * atomIndex, 4 * atomIndex + 3));
    }

    sum /= atomList.size();

    return sum;
  }

  /**
   * @brief Calculates the adjusted signed tetrahedron volume spanned by four
   *   sets of particle indices
   *
   * @param positions All positions
   * @param ligands Four particle index sets
   *
   * @note The volume is adjusted by V' = 6 * V
   *
   * @return The adjusted signed tetrahedron volume spanned by four sets of
   *   particle indices
   */
  inline static double adjustedSignedVolume(
    const Vector& positions,
    const ChiralityConstraint::SiteSequence& ligands
  ) {
    return (
      getAveragePos3D(positions, ligands[0])
      - getAveragePos3D(positions, ligands[3])
    ).dot(
      (
        getAveragePos3D(positions, ligands[1])
        - getAveragePos3D(positions, ligands[3])
      ).cross(
        getAveragePos3D(positions, ligands[2])
        - getAveragePos3D(positions, ligands[3])
      )
    );
  }

  /**
   * @brief Calculate the dihedral angle spanned by four particle index sets
   *
   * @param positions All positions
   * @param ligands The four particle index sets that span the dihedral
   *
   * @return The dihedral angle spanned by four particle index sets
   */
  inline static double dihedralAngle(
    const Vector& positions,
    const DihedralConstraint::SiteSequence& ligands
  ) {
    const dlib::vector<double, 3> alpha = getAveragePos3D(positions, ligands[0]);
    const dlib::vector<double, 3> beta = getAveragePos3D(positions, ligands[1]);
    const dlib::vector<double, 3> gamma = getAveragePos3D(positions, ligands[2]);
    const dlib::vector<double, 3> delta = getAveragePos3D(positions, ligands[3]);

    const dlib::vector<double, 3> g = gamma - beta;

    const dlib::vector<double, 3> m = (beta - alpha).cross(g);
    const dlib::vector<double, 3> n = g.cross(delta - gamma);

    return std::atan2(
      m.cross(n).dot(
        dlib::normalize(g)
      ),
      m.dot(n)
    );
  }
//!@}

  // State
  const unsigned N;
  const dlib::matrix<double, 0, 0>& squaredBounds;
  const std::vector<ChiralityConstraint>& chiralityConstraints;
  const std::vector<DihedralConstraint>& dihedralConstraints;
  bool compressFourthDimension = false;

  // Mutable state
  mutable double proportionChiralConstraintsCorrectSign = 0.0;
  mutable unsigned callCount = 0;

  // Constructor
  explicit ErrorFunctionValue(
    const dlib::matrix<double, 0, 0>& passSquaredBounds,
    const std::vector<ChiralityConstraint>& passChiralityConstraints,
    const std::vector<DihedralConstraint>& passDihedralConstraints
  ) : N(passSquaredBounds.nc()),
      squaredBounds(passSquaredBounds),
      chiralityConstraints(passChiralityConstraints),
      dihedralConstraints(passDihedralConstraints)
  {}

//!@name Information
//!@{
  /**
   * @brief Calculates the error function terms caused by distance errors
   * @param positions All positions
   * @return The sum of error function terms caused by distance error
   */
  double distanceError(const Vector& positions) const {
    double error = 0, lowerTerm, upperTerm;

    for(unsigned i = 0; i < N; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        const double upperBoundSquared = squaredBounds(i, j);
        const double lowerBoundSquared = squaredBounds(j, i);
        // Since i < j, upper Bound is (i, j), lower Bound is (j, i)
        assert(lowerBoundSquared <= upperBoundSquared);

        const double diffLength = dlib::length_squared(
          dlib::rowm(positions, dlib::range(4 * j, 4 * j + 3))
          - dlib::rowm(positions, dlib::range(4 * i, 4 * i + 3))
        );

        // Upper bound term
        upperTerm = diffLength / upperBoundSquared - 1;

        assert(!( // Ensure early continue logic is correct
          upperTerm > 0
          && (
            2 * lowerBoundSquared / (
              lowerBoundSquared + diffLength
            ) - 1
          ) > 0
        ));

        if(upperTerm > 0) {
          error += upperTerm * upperTerm;

          // If the upper term contributes, the lower certainly doesn't
          continue;
        }

        // Lower bound term
        lowerTerm = 2 * lowerBoundSquared / (
          lowerBoundSquared + diffLength
        ) - 1;

        if(lowerTerm > 0) {
          error += lowerTerm * lowerTerm;
        }
      }
    }

    return error;
  }

  /**
   * @brief Calculate the error function terms caused by chiral errors
   * @param positions All positions
   * @post Sets @p proportionChiralConstraintsCorrectSign
   * @return The sum of error function terms caused by chiral errors
   */
  double chiralError(const Vector& positions) const {
    unsigned nonZeroChiralityConstraints = 0;
    unsigned incorrectNonZeroChiralityConstraints = 0;

    double error = 0, volume, upperTerm, lowerTerm;

    for(const auto& constraint : chiralityConstraints) {
      volume = adjustedSignedVolume(positions, constraint.sites);

      if(std::fabs(constraint.lower + constraint.upper) > 1e-4) {
        ++nonZeroChiralityConstraints;
        if( // can this be simplified? -> sign bit XOR?
          ( volume < 0 && constraint.lower > 0)
          || (volume > 0 && constraint.lower < 0)
        ) {
          incorrectNonZeroChiralityConstraints += 1;
        }
      }

      // Upper bound term
      upperTerm = volume - constraint.upper;

      assert(!( // Ensure early continue logic is correct
        upperTerm > 0
        && constraint.lower - volume > 0
      ));

      if(upperTerm > 0) {
        error += upperTerm * upperTerm;

        // If the upper term contributes, the lower certainly doesn't
        continue;
      }

      // Lower bound term
      lowerTerm = constraint.lower - volume;
      if(lowerTerm > 0) {
        error += lowerTerm * lowerTerm;
      }
    }

    // Set proportionChiralConstraintsCorrectSign member
    if(nonZeroChiralityConstraints == 0) {
      proportionChiralConstraintsCorrectSign = 1;
    } else {
      proportionChiralConstraintsCorrectSign = static_cast<double>(
        nonZeroChiralityConstraints - incorrectNonZeroChiralityConstraints
      ) / nonZeroChiralityConstraints;
    }

    return error;
  }

  /**
   * @brief Calculate the error function term caused by dihedrals
   * @param positions All positions
   * @return The sum of error function terms caused by dihedrals
   */
  double dihedralError(const Vector& positions) const {
    double error = 0, phi, constraintSumHalved, term;

    for(const auto& constraint : dihedralConstraints) {
      phi = ErrorFunctionValue::dihedralAngle(positions, constraint.sites);

      constraintSumHalved = (constraint.lower + constraint.upper) / 2;

      if(phi < constraintSumHalved - M_PI) {
        phi += 2 * M_PI;
      } else if(phi > constraintSumHalved + M_PI) {
        phi -= 2 * M_PI;
      }

      term = std::fabs(phi - constraintSumHalved) - (constraint.upper - constraint.lower) / 2;

      if(term > 0) {
        error += term * term / 10;
      }
    }

    return error;
  }

  /**
   * @brief Calculates the error contribution due to the fourth spatial dimension
   * @note This should be calculated only if compress is true, where the
   *   refinement compresses out the fourth dimension
   * @param positions All positions
   * @return The sum of the squares of all fourth spatial dimension values
   */
  double extraDimensionError(const Vector& positions) const {
    double error = 0;
    for(unsigned i = 0; i < N; i++) {
      const double w = positions(4 * i + 3);
      error += w * w;
    }

    return error;
  }
//!@}

  /*!
   * @brief Function call operator for dlib optimization routine
   * Calculates all terms of the error function for a given set of positions
   */
  double operator() (const Vector& positions) const {
    ++callCount;

    assert(positions.size() % 4 == 0);
    assert(positions.size() / 4 == N);

    // Evaluation of addition proceeds left-to-right
    if(compressFourthDimension) {
      /* After chiral inversion, generally:
       * distance >> dihedral ~ chiral > extraDim
       */
      return (
        extraDimensionError(positions)
        + chiralError(positions)
        + distanceError(positions)
        + dihedralError(positions)
      );
    }

    /* Before chiral inversion, typically:
     * chiral >> dihedral ~ distance
     */
    return (
      distanceError(positions)
      + dihedralError(positions)
      + chiralError(positions)
    );
  }

  double calculateProportionChiralConstraintsCorrectSign(const Vector& positions) const {
    unsigned nonZeroChiralityConstraints = 0;
    unsigned incorrectNonZeroChiralityConstraints = 0;

    for(const auto& chiralityConstraint : chiralityConstraints) {
      /* Make sure that chirality constraints meant to cause coplanarity of
       * four indices aren't counted here - Their volume bounds are
       * usually so tight that they might be slightly outside in any given
       * refinement step. If they are outside their target values by a large
       * amount, that is caught by finalStructureAcceptable anyway.
       */
      if(std::fabs(chiralityConstraint.lower + chiralityConstraint.upper) > 1e-4) {
        nonZeroChiralityConstraints += 1;

        const double currentVolume = adjustedSignedVolume(
          positions,
          chiralityConstraint.sites
        );

        if( // can this be simplified? -> sign bit XOR?
          ( currentVolume < 0 && chiralityConstraint.lower > 0)
          || (currentVolume > 0 && chiralityConstraint.lower < 0)
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

  template<class Visitor>
  auto visitUnfulfilledConstraints(
    const DistanceBoundsMatrix& bounds,
    const Vector& positions,
    Visitor visitor
  ) const {
    // Check distance bounds deviations
    for(unsigned i = 0; i < bounds.N(); i++) {
      for(unsigned j = i + 1; j < bounds.N(); j++) {
        double ijDistance = dlib::length(
          getPos(positions, j) - getPos(positions, i)
        );

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
    for(const auto& constraint : chiralityConstraints) {
      double volume = adjustedSignedVolume(positions, constraint.sites);

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
      double phi = dihedralAngle(positions, constraint.sites);

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

/**
 * @brief Functor to calculate the gradient of the DG refinement error function
 *
 * Refinement progresses in two stages:
 * - Uncompressed, in which the spatial coordinates are free to expand into
 *   the fourth dimension. This is so that chiral constraints can invert
 *   more easily. Once all chiral constraints have the correct sign, we move
 *   on to the next stage.
 * - Compressed, in which the remaining distance bounds and chiral constraints
 *   are optimized further and the fourth dimension is "compressed" out.
 */
struct ErrorFunctionGradient {
public:
  // Typedefs
  using Vector = dlib::matrix<double, 0, 1>;

  // State
  const unsigned N;
  const dlib::matrix<double, 0, 0>& squaredBounds;
  const std::vector<ChiralityConstraint>& chiralityConstraints;
  const std::vector<DihedralConstraint>& dihedralConstraints;
  bool compressFourthDimension = false;

  // Constructor
  explicit ErrorFunctionGradient(
    const dlib::matrix<double, 0, 0>& passSquaredBounds,
    const std::vector<ChiralityConstraint>& passChiralityConstraints,
    const std::vector<DihedralConstraint>& passDihedralConstraints
  ) : N(passSquaredBounds.nc()),
      squaredBounds(passSquaredBounds),
      chiralityConstraints(passChiralityConstraints),
      dihedralConstraints(passDihedralConstraints)
  {}

//!@name Information
//!@{
  /**
   * @brief A reference, unoptimized implementation of the A term (The first
   * distance bounds error term)
   *
   * @param positions All positions
   *
   * @return A vector containing only the contributions of the A term to the
   *   gradient
   */
  Vector referenceA(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    for(unsigned alpha = 0; alpha < N; ++alpha) {
      for(unsigned i = 0; i < N; ++i) {
        if(i == alpha) { // skip identical indices
          continue;
        }

        const Vector diff = ErrorFunctionValue::getPos(positions, alpha) - ErrorFunctionValue::getPos(positions, i);

        const double upperBoundSquared = squaredBounds(std::min(i, alpha), std::max(i, alpha));

        const double upperTerm = dlib::length_squared(diff) / upperBoundSquared - 1;

        if(upperTerm > 0) {
          dlib::set_rowm(
            gradient,
            dlib::range(4 * alpha, 4 * alpha + 3)
          ) += 4 * diff * upperTerm / upperBoundSquared;
        }
      }
    }

    return gradient;
  }

  /**
   * @brief A reference, unoptimized implementation of the B term (The second
   * distance bounds error term)
   *
   * @param positions All positions
   *
   * @return A vector containing only the contributions of the B term to the
   *   gradient
   */
  Vector referenceB(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    for(unsigned alpha = 0; alpha < N; ++alpha) {
      for(unsigned i = 0; i < N; ++i) {
        if(i == alpha) { // skip identical indices
          continue;
        }

        const Vector diff = ErrorFunctionValue::getPos(positions, i) - ErrorFunctionValue::getPos(positions, alpha);

        const double lowerBoundSquared = squaredBounds(
          std::max(i, alpha),
          std::min(i, alpha)
        );

        const double quotient = lowerBoundSquared + dlib::length_squared(diff);

        const double lowerTerm = 2 * lowerBoundSquared / quotient - 1;

        if(lowerTerm > 0) {
          dlib::set_rowm(
            gradient,
            dlib::range(4 * alpha, 4 * alpha + 3)
          ) += 8 * lowerBoundSquared * diff * lowerTerm / (
            quotient * quotient
          );
        }
      }
    }

    return gradient;
  }

  /**
   * @brief Adds a chiral errors contribution to the gradient vector
   *
   * @param gradient A reference to the gradient vector
   * @param positions All positions
   * @param constraint A particular chiral constraint
   * @param factor A prefactor with which to multiply contributions
   */
  void referenceAddChiralContribution(
    Vector& gradient,
    const Vector& positions,
    const ChiralityConstraint& constraint,
    const double factor = 1
  ) const {
    // Precalculate repeated expressions
    const dlib::vector<double, 3> alphaMinusDelta = (
      ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[0])
      - ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[3])
    );

    const dlib::vector<double, 3> betaMinusDelta = (
      ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[1])
      - ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[3])
    );

    const dlib::vector<double, 3> gammaMinusDelta = (
      ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[2])
      - ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[3])
    );

    // Specific to deltaI only but still repeated
    const dlib::vector<double, 3> alphaMinusGamma = (
      ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[0])
      - ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[2])
    );

    const dlib::vector<double, 3> betaMinusGamma = (
      ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[1])
      - ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[2])
    );

    const dlib::vector<double, 3> iContribution =
      factor * betaMinusDelta.cross(gammaMinusDelta)
      / constraint.sites[0].size();

    const dlib::vector<double, 3> jContribution =
      factor * gammaMinusDelta.cross(alphaMinusDelta)
      / constraint.sites[1].size();

    const dlib::vector<double, 3> kContribution =
      factor * alphaMinusDelta.cross(betaMinusDelta)
      / constraint.sites[2].size();

    const dlib::vector<double, 3> lContribution =
      factor * betaMinusGamma.cross(alphaMinusGamma)
      / constraint.sites[3].size();

    for(const auto alphaI : constraint.sites[0]) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * alphaI, 4 * alphaI + 2)
      ) += iContribution;
    }

    for(const auto betaI : constraint.sites[1]) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * betaI, 4 * betaI + 2)
      ) += jContribution;
    }

    for(const auto gammaI : constraint.sites[2]) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * gammaI, 4 * gammaI + 2)
      ) += kContribution;
    }

    for(const auto deltaI : constraint.sites[3]) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * deltaI, 4 * deltaI + 2)
      ) += lContribution;
    }
  }

  /**
   * @brief A reference, unoptimized implementation of the C term (all chiral
   * error terms)
   *
   * @param positions All positions
   *
   * @return A vector containing only the contributions of the C term to the
   *   gradient
   */
  Vector referenceC(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    for(const auto& constraint : chiralityConstraints) {
      const double volume = ErrorFunctionValue::adjustedSignedVolume(positions, constraint.sites);
      const double upperTerm = volume - constraint.upper;
      const double lowerTerm = constraint.lower - volume;

      if(upperTerm > 0 || lowerTerm > 0) {
        referenceAddChiralContribution(
          gradient,
          positions,
          constraint,
          2 * (std::max(0.0, upperTerm) - std::max(0.0, lowerTerm))
        );
      }
    }

    return gradient;
  }

  /**
   * @brief A reference, unoptimized implementation of the D term (compression)
   *
   * @param positions All positions
   *
   * @return A vector containing only the contributions of the D term to the
   *   gradient
   */
  Vector referenceD(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    if(compressFourthDimension) {
      for(unsigned i = 0; i < N; ++i) {
        gradient(4 * i + 3) += 2 * positions(4 * i + 3);
      }
    }

    return gradient;
  }

  Vector referenceDihedral(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    // Currently no slower reference implementation of the same gradient

    addDihedralContributions(gradient, positions);

    return gradient;
  }

  /**
   * @brief A reference, unoptimized implementation of the gradient calculation
   *
   * Exists solely to double-check the optimized calculations.
   */
  Vector reference(const Vector& positions) const {
    return (
      referenceA(positions)
      + referenceB(positions)
      + referenceC(positions)
      + referenceD(positions)
      + referenceDihedral(positions)
    );
  }

  /**
   * @brief Adds the contribution of a pair of atoms to the gradient
   *
   * @param positions All positions
   * @param gradient A reference to the gradient under construction
   * @param lowerBoundSquared The square of the lower bound
   * @param upperBoundSquared The square of the upper bound
   * @param alpha
   * @param i
   */
  inline void gradientDistanceContribution(
    const Vector& positions,
    Vector& gradient,
    const double lowerBoundSquared,
    const double upperBoundSquared,
    const unsigned alpha,
    const unsigned i
  ) const {
    assert(lowerBoundSquared <= upperBoundSquared);

    // For both
    const Vector alphaMinusI = (
      dlib::rowm(positions, dlib::range(4 * alpha, 4 * alpha + 3))
      - dlib::rowm(positions, dlib::range(4 * i, 4 * i + 3))
    );

    const double diffLength = dlib::length_squared(alphaMinusI);

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
      const Vector f = 4 * alphaMinusI * upperTerm / upperBoundSquared;

      dlib::set_rowm(
        gradient,
        dlib::range(4 * alpha, 4 * alpha + 3)
      ) += f;

      dlib::set_rowm(
        gradient,
        dlib::range(4 * i, 4 * i + 3)
      ) -= f;

      /* if the upper term is set, there is no way the lower term needs to be
       * applied -> early return
       */
      return;
    }

    // Lower term
    const double quotient = lowerBoundSquared + diffLength;
    const double lowerTerm = 2 * lowerBoundSquared / quotient - 1;

    if(lowerTerm > 0) {
      const Vector g = 8 * lowerBoundSquared * alphaMinusI * lowerTerm / (
        quotient * quotient
      );

      /* We use -= because the lower term needs the position vector
       * difference (i - alpha), so we reuse alphaMinusI and just subtract
       * from the gradient instead of adding to it
       */
      dlib::set_rowm(
        gradient,
        dlib::range(4 * alpha, 4 * alpha + 3)
      ) -= g;

      dlib::set_rowm(
        gradient,
        dlib::range(4 * i, 4 * i + 3)
      ) += g;
    }
  }

  /**
   * @brief Adds the contribution of all dihedral constraints to a gradient vector
   *
   * @param gradient The gradient vector to add dihedral constraint
   *   contributions to
   * @param positions The current positions vector
   */
  inline void addDihedralContributions(
    Vector& gradient,
    const Vector& positions
  ) const {
    for(const DihedralConstraint& constraint : dihedralConstraints) {
      const dlib::vector<double, 3> alpha = ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[0]);
      const dlib::vector<double, 3> beta = ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[1]);
      const dlib::vector<double, 3> gamma = ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[2]);
      const dlib::vector<double, 3> delta = ErrorFunctionValue::getAveragePos3D(positions, constraint.sites[3]);

      const dlib::vector<double, 3> f = alpha - beta;
      const dlib::vector<double, 3> g = beta - gamma;
      const dlib::vector<double, 3> h = delta - gamma;

      dlib::vector<double, 3> a = f.cross(g);
      dlib::vector<double, 3> b = h.cross(g);

      // a and b can be zero-length-vectors!

      // All of the permutations yield identical expressions within atan2 aside from -g
      double phi = std::atan2(
        a.cross(b).dot(
          dlib::normalize(-g)
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

      // Multiply with sgn (w)
      h_phi *= static_cast<int>(0.0 < w_phi) - static_cast<int>(w_phi < 0.0);

      // Multiply in the factors 2 * 1/10 -> 1/5
      h_phi /= 5;

      // Precompute some reused expressions
      const double gLength = dlib::length(g);
      if(gLength == 0) {
        // All contributions would be NaNs
        continue;
      }

      const double aLengthSq = dlib::length_squared(a);
      if(aLengthSq > 0) {
        a /= aLengthSq;
      }
      const double bLengthSq = dlib::length_squared(b);
      if(bLengthSq > 0) {
        b /= bLengthSq;
      }
      const double fDotG = f.dot(g);
      const double gDotH = g.dot(h);

      const dlib::vector<double, 3> iContribution = -(h_phi / constraint.sites[0].size()) * gLength * a;

      const dlib::vector<double, 3> jContribution = (h_phi / constraint.sites[1].size()) * (
        (gLength + fDotG / gLength) * a
        - (gDotH / gLength) * b
      );

      const dlib::vector<double, 3> kContribution = (h_phi / constraint.sites[2].size()) * (
        (gDotH / gLength - gLength) * b
        - (fDotG / gLength) * a
      );

      const dlib::vector<double, 3> lContribution = (h_phi / constraint.sites[3].size()) * gLength * b;

      if(
        !dlib::is_finite(iContribution)
        || !dlib::is_finite(jContribution)
        || !dlib::is_finite(kContribution)
        || !dlib::is_finite(lContribution)
      ) {
        throw std::out_of_range("Encountered non-finite dihedral contributions");
      }

      for(const AtomIndex alphaConstitutingIndex : constraint.sites[0]) {
        dlib::set_rowm(
          gradient,
          dlib::range(4 * alphaConstitutingIndex, 4 * alphaConstitutingIndex + 2)
        ) += iContribution;
      }

      for(const AtomIndex betaConstitutingIndex: constraint.sites[1]) {
        dlib::set_rowm(
          gradient,
          dlib::range(4 * betaConstitutingIndex, 4 * betaConstitutingIndex + 2)
        ) += jContribution;
      }

      for(const AtomIndex gammaConstitutingIndex : constraint.sites[2]) {
        dlib::set_rowm(
          gradient,
          dlib::range(4 * gammaConstitutingIndex, 4 * gammaConstitutingIndex + 2)
        ) += kContribution;
      }

      for(const AtomIndex deltaConstitutingIndex : constraint.sites[3]) {
        dlib::set_rowm(
          gradient,
          dlib::range(4 * deltaConstitutingIndex, 4 * deltaConstitutingIndex + 2)
        ) += lContribution;
      }
    }
  }

//!@}

  /*!
   * @brief Function call operator for dlib optimization routine
   * Calculates all terms of the error function gradient for a given set of positions
   */
  Vector operator() (const Vector& positions) const {
    assert(positions.size() % 4 == 0);
    assert(positions.size() / 4 == N);

    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    // Distance gradient contributions (A and B)
    for(unsigned i = 0; i < N - 1; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        gradientDistanceContribution(
          positions,
          gradient,
          squaredBounds(j, i), // lower bound
          squaredBounds(i, j), // upper bound
          i,
          j
        );
      }
    }

    // Chirality gradient contributions (C)
    for(const auto& constraint : chiralityConstraints) {
      const double volume = ErrorFunctionValue::adjustedSignedVolume(positions, constraint.sites);

      // It may not occur that both terms contribute!
      assert(!(volume - constraint.upper > 0 && constraint.lower - volume > 0));

      const double factor = 2 * (
        std::max(0.0, volume - constraint.upper)
        - std::max(0.0, constraint.lower - volume)
      );

      /* Make sure computing all cross products is worth it by checking that
       * one of both terms is actually greater than 0
       */
      if(factor != 0.0) {
        referenceAddChiralContribution(
          gradient,
          positions,
          constraint,
          factor
        );
      }
    }

    // Add dihedral errors
    addDihedralContributions(gradient, positions);

    // Fourth dimension contribution (E)
    if(compressFourthDimension) {
      for(unsigned i = 0; i < N; ++i) {
        gradient(4 * i + 3) += 2 * positions(4 * i + 3);
      }
    }

    return gradient;
  }
};

} // namespace DistanceGeometry
} // namespace molassembler
} // namespace Scine

#endif
