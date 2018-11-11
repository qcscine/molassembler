// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_DG_DLIB_REFINEMENT_PROBLEM_H
#define INCLUDE_MOLASSEMBLER_DG_DLIB_REFINEMENT_PROBLEM_H

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"

#include <dlib/optimization.h>

/*! @file
 *
 * @brief Specifies the refinement minimization problem for the dlib library.
 */

namespace molassembler {

namespace DistanceGeometry {

namespace errfDetail {

using Vector = dlib::matrix<double, 0, 1>;

/**
 * @brief Fetches a three-dimensional position sub-part of the positions vector
 *
 * @param positions All positions
 * @param i The index of the particle whose position to extract
 *
 * @return The position of a particle in the positions
 */
inline dlib::vector<double, 3> getPos3D(const Vector& positions, const unsigned i) {
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
inline Vector getPos(const Vector& positions, const unsigned i) {
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
inline dlib::vector<double, 3> getAveragePos3D(
  const Vector& positions,
  const ChiralityConstraint::AtomListType& atomList
) {
  if(atomList.size() == 1) {
    return getPos3D(positions, atomList.front());
  }

  dlib::vector<double, 3> sum {0.0, 0.0, 0.0};

  for(const auto& atomIndex : atomList) {
    sum += getPos3D(positions, atomIndex);
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
inline Vector getAveragePos(
  const Vector& positions,
  const ChiralityConstraint::AtomListType& atomList
) {
  if(atomList.size() == 1) {
    return getPos(positions, atomList.front());
  }

  Vector sum = dlib::zeros_matrix<double>(4, 1);

  for(const auto& atomIndex : atomList) {
    sum += getPos(positions, atomIndex);
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
inline double adjustedSignedVolume(
  const Vector& positions,
  const ChiralityConstraint::LigandSequence& ligands
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
inline double dihedralAngle(
  const Vector& positions,
  const DihedralConstraint::LigandSequence& ligands
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

/**
 * @brief Calculates the proportion of chirality constraints with the correct
 *   sign
 *
 * @param chiralityConstraints The list of chirality constraints
 * @param positions All positions
 *
 * @return The ratio of chirality constraints with the correct sign to the
 *   total number of chirality constraints
 */
double proportionChiralityConstraintsCorrectSign(
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  const Vector& positions
);

/**
 * @brief Decides whether the final structure from a refinement is acceptable
 *
 * A final structure is acceptable if
 * - All distance bounds are within 0.5 of either the lower or upper boundary
 * - All chiral constraints are within 0.5 of either the lower or upper boundary
 *
 * @param bounds The distance bounds
 * @param chiralityConstraints All chirality constraints
 * @param positions The final positions from a refinement
 *
 * @return Whether the final structure is acceptable
 */
bool finalStructureAcceptable(
  const DistanceBoundsMatrix& bounds,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  const Vector& positions
);

/**
 * @brief Writes messages to the log explaining why a structure was rejected
 *
 * @see finalStructureAcceptable for criteria on acceptability
 *
 * A final structure is acceptable if
 * - All distance bounds are within 0.5 of either the lower or upper boundary
 * - All chiral constraints are within 0.5 of either the lower or upper boundary
 *
 * @param bounds The distance bounds
 * @param chiralityConstraints All chirality constraints
 * @param positions The final positions from a refinement
 */
void explainAcceptanceFailure(
  const DistanceBoundsMatrix& bounds,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  const Vector& positions
);

/**
 * @brief Writes messages to the log explaining the final contributions to
 *   the error function
 *
 * @param bounds The distance bounds
 * @param chiralityConstraints All chirality constraints
 * @param positions The final positions from a refinement
 */
void explainFinalContributions(
  const DistanceBoundsMatrix& bounds,
  const std::vector<ChiralityConstraint>& chiralityConstraints,
  const Vector& positions
);

} // namespace errfDetail

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
 *
 * @tparam compress Whether to compress out the fourth spatial dimension (i.e.
 *   penalize non-zero values of the fourth dimension).
 */
template<bool compress>
struct errfValue {
  // Typedefs
  using Vector = dlib::matrix<double, 0, 1>;

  // State
  const unsigned N;
  const dlib::matrix<double, 0, 0>& squaredBounds;
  const std::vector<ChiralityConstraint>& chiralityConstraints;
  const std::vector<DihedralConstraint>& dihedralConstraints;

  // Constructor
  explicit errfValue(
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
          errfDetail::getPos(positions, j) - errfDetail::getPos(positions, i)
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
   * @return The sum of error function terms caused by chiral errors
   */
  double chiralError(const Vector& positions) const {
    double error = 0, volume, upperTerm, lowerTerm;

    for(const auto& constraint : chiralityConstraints) {
      volume = errfDetail::adjustedSignedVolume(positions, constraint.sites);

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
      phi = dihedralAngle(positions, constraint.sites);

      constraintSumHalved = (constraint.lower + constraint.upper) / 2;

      if(phi < constraintSumHalved - M_PI) {
        phi += 2 * M_PI;
      } else if(phi > constraintSumHalved + M_PI) {
        phi -= 2 * M_PI;
      }

      term = std::fabs(phi - constraintSumHalved) - (constraint.upper - constraint.lower) / 2;

      if(term > 0) {
        error += term * term;
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
    assert(positions.size() % 4 == 0);
    assert(positions.size() / 4 == N);

    // Evaluation of addition proceeds left-to-right
    if(compress) {
      // After chiral inversion, generally: distance >> chiral > extraDim
      return (
        extraDimensionError(positions)
        + chiralError(positions)
        + distanceError(positions)
      );
    }

    // Before chiral inversion, typically: chiral >> distance
    return (
      distanceError(positions)
      + chiralError(positions)
    );
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
 *
 * @tparam compress Whether to compress out the fourth spatial dimension (i.e.
 *   penalize non-zero values of the fourth dimension).
 */
template<bool compress>
struct errfGradient {
public:
  // Typedefs
  using Vector = dlib::matrix<double, 0, 1>;

  // State
  const unsigned N;
  const dlib::matrix<double, 0, 0>& squaredBounds;
  const std::vector<ChiralityConstraint>& chiralityConstraints;
  const std::vector<DihedralConstraint>& dihedralConstraints;

  // Constructor
  explicit errfGradient(
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

    for(unsigned alpha = 0; alpha < N; alpha++) {
      for(unsigned i = 0; i < N; i++) {
        if(i == alpha) { // skip identical indices
          continue;
        }

        const Vector diff = errfDetail::getPos(positions, alpha) - errfDetail::getPos(positions, i);

        const double upperTerm = dlib::length_squared(diff) / squaredBounds(
          std::min(i, alpha),
          std::max(i, alpha)
        ) - 1;

        if(upperTerm > 0) {
          dlib::set_rowm(
            gradient,
            dlib::range(4 * alpha, 4 * alpha + 3)
          ) += 4 * diff * upperTerm / squaredBounds(
            std::min(i, alpha),
            std::max(i, alpha)
          );
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

    for(unsigned alpha = 0; alpha < N; alpha++) {
      for(unsigned i = 0; i < N; i++) {
        if(i == alpha) { // skip identical indices
          continue;
        }

        const Vector diff = errfDetail::getPos(positions, i) - errfDetail::getPos(positions, alpha);

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
      errfDetail::getAveragePos3D(positions, constraint.sites[0])
      - errfDetail::getAveragePos3D(positions, constraint.sites[3])
    );

    const dlib::vector<double, 3> betaMinusDelta = (
      errfDetail::getAveragePos3D(positions, constraint.sites[1])
      - errfDetail::getAveragePos3D(positions, constraint.sites[3])
    );

    const dlib::vector<double, 3> gammaMinusDelta = (
      errfDetail::getAveragePos3D(positions, constraint.sites[2])
      - errfDetail::getAveragePos3D(positions, constraint.sites[3])
    );

    // Specific to deltaI only but still repeated
    const dlib::vector<double, 3> alphaMinusGamma = (
      errfDetail::getAveragePos3D(positions, constraint.sites[0])
      - errfDetail::getAveragePos3D(positions, constraint.sites[2])
    );

    const dlib::vector<double, 3> betaMinusGamma = (
      errfDetail::getAveragePos3D(positions, constraint.sites[1])
      - errfDetail::getAveragePos3D(positions, constraint.sites[2])
    );

    for(const auto alphaI : constraint.sites[0]) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * alphaI, 4 * alphaI + 2)
      ) += factor * betaMinusDelta.cross(gammaMinusDelta) / constraint.sites[0].size();
    }

    for(const auto betaI : constraint.sites[1]) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * betaI, 4 * betaI + 2)
      ) += factor * gammaMinusDelta.cross(alphaMinusDelta) / constraint.sites[1].size();
    }

    for(const auto gammaI : constraint.sites[2]) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * gammaI, 4 * gammaI + 2)
      ) += factor * alphaMinusDelta.cross(betaMinusDelta) / constraint.sites[2].size();
    }

    for(const auto deltaI : constraint.sites[3]) {
      dlib::set_rowm(
        gradient,
        dlib::range(4 * deltaI, 4 * deltaI + 2)
      ) += factor * betaMinusGamma.cross(alphaMinusGamma) / constraint.sites[3].size();
    }
  }

  /**
   * @brief A reference, unoptimized implementation of the C term (The first
   * chiral bounds error term)
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
      const double volume = errfDetail::adjustedSignedVolume(positions, constraint.sites);
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
   * @brief A reference, unoptimized implementation of the D term (The second
   * chiral bounds error term)
   *
   * @param positions All positions
   *
   * @return A vector containing only the contributions of the D term to the
   *   gradient
   */
  Vector referenceD(const Vector& positions) const {
    Vector gradient(positions.size());
    dlib::set_all_elements(gradient, 0);

    if(compress) {
      for(unsigned i = 0; i < N; i++) {
        gradient(4 * i + 3) += 2 * positions(4 * i + 3);
      }
    }

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
      errfDetail::getPos(positions, alpha)
      - errfDetail::getPos(positions, i)
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
      dlib::set_rowm(
        gradient,
        dlib::range(4 * alpha, 4 * alpha + 3)
      ) += 4 * alphaMinusI * upperTerm / upperBoundSquared;

      /* if the upper term is set, there is no way the lower term needs to be
       * applied -> early return
       */
      return;
    }

    // Lower term
    const double quotient = lowerBoundSquared + diffLength;
    const double lowerTerm = 2 * lowerBoundSquared / quotient - 1;

    if(lowerTerm > 0) {
      /* We use -= because the lower term needs the position vector
       * difference (i - alpha), so we reuse alphaMinusI and just subtract
       * from the gradient instead of adding to it
       */
      dlib::set_rowm(
        gradient,
        dlib::range(4 * alpha, 4 * alpha + 3)
      ) -= 8 * lowerBoundSquared * alphaMinusI * lowerTerm / (
        quotient * quotient
      );
    }
  }

  inline void addDihedralConstributions(
    const Vector& positions,
    Vector& gradient
  ) const {
    double phi, constraintSumHalved, w_phi, h_phi;

    for(const auto& constraint : dihedralConstraints) {
      const dlib::vector<double, 3> alpha = getAveragePos3D(positions, constraint.sites[0]);
      const dlib::vector<double, 3> beta = getAveragePos3D(positions, constraint.sites[1]);
      const dlib::vector<double, 3> gamma = getAveragePos3D(positions, constraint.sites[2]);
      const dlib::vector<double, 3> delta = getAveragePos3D(positions, constraint.sites[3]);

      const dlib::vector<double, 3> f = beta - alpha;
      const dlib::vector<double, 3> g = gamma - beta;
      const dlib::vector<double, 3> h = delta - gamma;

      const dlib::vector<double, 3> a = f.cross(g);
      const dlib::vector<double, 3> b = g.cross(h);

      phi = std::atan2(
        a.cross(b).dot(
          dlib::normalize(g)
        ),
        a.dot(b)
      );

      constraintSumHalved = (constraint.lower + constraint.upper) / 2;

      if(phi < constraintSumHalved - M_PI) {
        phi += 2 * M_PI;
      } else if(phi > constraintSumHalved + M_PI) {
        phi -= 2 * M_PI;
      }

      w_phi = phi - constraintSumHalved;

      h_phi = std::fabs(w_phi) - (constraint.upper - constraint.lower) / 2;

      /* "Apply the max function": If h <= 0, then the max function yields zero,
       * so we can skip this
       */
      if(h_phi <= 0) {
        continue;
      }

      // Multiply with sgn (w) -> (-1) if
      if(w_phi < 0) {
        h_phi *= -1;
      }

      // Add in the prefactor
      h_phi *= 2.0;

      // Precompute some reused expressions
      const double gNorm = dlib::norm(g);
      const dlib::vector<double, 3> aAdjusted = a / dlib::length_squared(a);
      const dlib::vector<double, 3> bAdjusted = b / dlib::length_squared(b);
      std::array<unsigned, 4> siteSizes;
      for(unsigned i = 0; i < 4; ++i) {
        siteSizes[i] = constraint.sites[i].size();
      }
      const double fDotG = f.dot(g);
      const double gDotH = g.dot(h);

      for(const auto& alphaConstitutingIndex : constraint.sites[0]) {
        dlib::set_rowm(
          gradient,
          dlib::range(4 * alphaConstitutingIndex, 4 * alphaConstitutingIndex + 2)
        ) -= gNorm * aAdjusted / siteSizes[0];
      }

      for(const auto& betaConstitutingIndex: constraint.sites[1]) {
        dlib::set_rowm(
          gradient,
          dlib::range(4 * betaConstitutingIndex, 4 * betaConstitutingIndex + 2)
        ) += (
          gNorm * aAdjusted
          + fDotG * aAdjusted / gNorm
          - gDotH * bAdjusted / gNorm
        ) / siteSizes[1];
      }

      for(const auto& gammaConstitutingIndex : constraint.sites[3]) {
        dlib::set_rowm(
          gradient,
          dlib::range(4 * gammaConstitutingIndex, 4 * gammaConstitutingIndex + 2)
        ) += (
          - gNorm * bAdjusted
          + gDotH * bAdjusted / gNorm
          - fDotG * aAdjusted / gNorm
        ) / siteSizes[2];
      }

      for(const auto& deltaConstitutingIndex : constraint.sites[3]) {
        dlib::set_rowm(
          gradient,
          dlib::range(4 * deltaConstitutingIndex, 4 * deltaConstitutingIndex + 2)
        ) += gNorm * bAdjusted / siteSizes[3];
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
    for(unsigned alpha = 0; alpha < N; alpha++) {
      for(unsigned i = 0; i < alpha; i++) {
        // Here, i < alpha, so lower bound is (alpha, i), upper bound is (i, alpha)
        gradientDistanceContribution(
          positions,
          gradient,
          squaredBounds(alpha, i),
          squaredBounds(i, alpha),
          alpha,
          i
        );
      }

      for(unsigned i = alpha + 1; i < N; i++) {
        // Here, i > alpha, so lower bound is (i, alpha), upper bound is (alpha, i)
        gradientDistanceContribution(
          positions,
          gradient,
          squaredBounds(i, alpha),
          squaredBounds(alpha, i),
          alpha,
          i
        );
      }
    }

    // Chirality gradient contributions (C)
    for(const auto& constraint : chiralityConstraints) {
      const double volume = errfDetail::adjustedSignedVolume(positions, constraint.sites);

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

    // Fourth dimension contribution (E)
    if(compress) {
      for(unsigned i = 0; i < N; i++) {
        gradient(4 * i + 3) += 2 * positions(4 * i + 3);
      }
    }

    return gradient;
  }
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
