#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H

#include "chemical_symmetries/Names.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/TinySet.h"

namespace stereopermutation {

namespace detail {

template<typename T>
bool nextCombinationPermutation(
  std::vector<T>& toPermute,
  const std::vector<T>& limits
) {
  assert(toPermute.size() == limits.size());
  const unsigned cols = toPermute.size();

  // Check if all columns are full
  bool allFull = true;
  for(unsigned i = 0; i < cols; ++i) {
    if(toPermute[i] != limits[i]) {
      allFull = false;
      break;
    }
  }

  if(allFull) {
    return false;
  } else {
    // Make next permutation
    for(int i = cols - 1; i >= 0; --i) {
      if(toPermute[i] == limits[i]) {
        toPermute[i] = 0;
      } else {
        ++toPermute[i];
        return true;
      }
    }

    return true;
  }
}

template<typename T>
std::pair<T, T> makeOrderedPair(T a, T b) {
  std::pair<T, T> pair {
    std::move(a),
    std::move(b)
  };

  if(pair.second < pair.first) {
    std::swap(pair.first, pair.second);
  }

  return pair;
}

} // namespace detail

class Composite {
public:
  enum class Orientation {
    Ecliptic,
    Staggered
  };

  using DihedralTuple = std::tuple<unsigned, unsigned, double>;

private:
  Symmetry::Name _left, _right;
  unsigned _leftFusedPosition, _rightFusedPosition;

  std::vector<
    std::vector<DihedralTuple>
  > _stereopermutations;

  /*!
   * Calculates the angle between two substituents that have the same angle from
   * the bound symmetry position in a plane perpendicular to their shared axis.
   *
   * The axis is defined by the bound symmetry position and the central point of
   * the symmetry. Symmetry positions with identical angles from the bound
   * symmetry position are situated together in a plane perpendicular to the
   * axis. The angle between these symmetry positions via the intersection with
   * the axis in that plane is calculated by this function.
   *
   * \warning Do not call this with angleFromBoundSymmetryPosition == M_PI as
   * the geometrical idea collapses.
   */
  static double perpendicularSubstituentAngle(
    const double angleFromBoundSymmetryPosition,
    const double angleBetweenSubstituents
  );

public:
  //! 1° absolute tolerance floating point comparison helper for angle groups
  static constexpr temple::floating::ExpandedAbsoluteEqualityComparator<double> fpComparator {
    temple::Math::toRadians(1.0)
  };

  //! Returns the
  static std::vector<unsigned> identityRotation(const Symmetry::Name symmetry) ;

  /*! Generates a symmetry rotation subject to constraints
   *
   * \param symmetryName The symmetry in which the rotation is sought
   * \param fixedSymmetryPosition The symmetry position in symmetryName that is
   *   to be kept constant
   * \param changedPositions The symmetry positions which must change in the
   *   sought rotation
   */
  static std::vector<unsigned> generateRotation(
    const Symmetry::Name symmetryName,
    const unsigned fixedSymmetryPosition,
    const std::vector<unsigned>& changedPositions
  );

  static std::vector<unsigned> rotation(
    const Symmetry::Name symmetryName,
    const unsigned fixedSymmetryPosition,
    const std::vector<unsigned>& perpendicularPlanePositions
  );

  struct AngleGroup {
    double angle;
    std::vector<unsigned> symmetryPositions;
  };

  //! Collects all coplanar indices that are closest to the fused symmetry position
  static AngleGroup smallestAngleGroup(
    const Symmetry::Name symmetryName,
    const unsigned fusedSymmetryPosition
  );

  using PerpendicularAngleGroups = std::vector<
    std::pair<
      temple::TinyUnorderedSet<double>,
      std::vector<
        std::pair<unsigned, unsigned>
      >
    >
  >;

  //! Creates sets of within-group cross angles in the perpendicular plane
  static PerpendicularAngleGroups inGroupAngles(
    const AngleGroup& angleGroup,
    const Symmetry::Name symmetryName
  );

  Composite() = default;

  Composite(
    const Symmetry::Name left,
    const Symmetry::Name right,
    const unsigned leftFusedPosition,
    const unsigned rightFusedPosition,
    const std::vector<char>& leftCharacters,
    const std::vector<char>& rightCharacters
  );

  //! Returns the number of permutations for this Composite
  unsigned permutations() const;

  /*! Returns a set of dihedrals for a particular permutation
   *
   * \note The first two elements of each tuple specify the symmetry position
   * within that side's symmetry. The first element is for the left symmetry,
   * the second for the right symmetry.
   */
  const std::vector<DihedralTuple>& dihedrals(unsigned permutationIndex) const;

  inline Symmetry::Name left() const {
    return _left;
  }

  inline Symmetry::Name right() const {
    return _right;
  }

  inline unsigned leftFusedPosition() const {
    return _leftFusedPosition;
  }

  inline unsigned rightFusedPosition() const {
    return _rightFusedPosition;
  }

  bool operator == (const Composite& other) const;
  bool operator != (const Composite& other) const;
};

} // namespace stereopermutation

#endif
