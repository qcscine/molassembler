#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H

#include "chemical_symmetries/Names.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/TinySet.h"
#include "temple/OrderedPair.h"

namespace stereopermutation {

// TODO these are unneeded here, just move to impl file
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
  struct AngleGroup {
    double angle;
    std::vector<unsigned> symmetryPositions;
    bool isotropic;
  };

  struct OrientationState {
    //! The symmetry of either positional symmetry
    Symmetry::Name symmetry;
    //! The symmetry position of the local symmetry that the other is fused at
    unsigned fusedPosition;
    //! Abstract ranking-characters of the ligands at their symmetry positions
    std::vector<char> characters;
    /*! An identifier to the symmetry source
     *
     * Since OrientationState is used internally in an ordered pair which may
     * swap elements on changes, an identifier can be useful in reassociating
     * the OrientationState with its data source.
     */
    std::size_t identifier;

    OrientationState(
      Symmetry::Name passSymmetry,
      unsigned passFusedPosition,
      std::vector<char> passCharacters,
      std::size_t passIdentifier
    );

    void applyCharacterRotation(const std::vector<unsigned>& rotation);

    //! Smallest symmetry position from the same group as the fused position
    unsigned lowestEqualPositionInSymmetry() const;

    //! Calculates the required reduction mapping to the canonical form
    std::vector<unsigned> findReductionMapping(unsigned reducedFusedPosition) const;

    /*! Transforms the OrientationState to a canonical form
     *
     * Transforms the OrientationState by applying a reduction mapping to the
     * smallest symmetry position from the same group as the fused position.
     * Returns the mapping needed to revert the OrientationState back to its
     * original data values.
     */
    [[nodiscard]] std::vector<unsigned> transformToCanonical();

    //! Reverts the OrientationState to non-canonical form
    void revert(const std::vector<unsigned>& reversionMapping);

    //! Collects all coplanar indices that are closest to the fused symmetry position
    AngleGroup smallestAngleGroup() const;

    //! Full member lexicographical comparison in order of declaration
    bool operator < (const OrientationState& other) const;
    //! Full member lexicographical comparison in order of declaration
    bool operator == (const OrientationState& other) const;
  };

  using DihedralTuple = std::tuple<unsigned, unsigned, double>;
  using PermutationsList = std::vector<
    std::vector<DihedralTuple>
  >;

private:
  temple::OrderedPair<OrientationState> _orientations;

  PermutationsList _stereopermutations;

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

  Composite(OrientationState first, OrientationState second);

  //! Returns the number of permutations for this Composite
  unsigned permutations() const;

  /*! Returns a set of dihedrals for a particular permutation
   *
   * \note The first two elements of each tuple specify the symmetry position
   * within that side's symmetry. The first element is for the left symmetry,
   * the second for the right symmetry.
   */
  const std::vector<DihedralTuple>& dihedrals(unsigned permutationIndex) const;

  //! Returns the orientation state of the composite
  const temple::OrderedPair<OrientationState>& orientations() const;

  //! Through-iteration of the generated dihedral permutations
  PermutationsList::const_iterator begin() const;
  PermutationsList::const_iterator end() const;

  bool operator == (const Composite& other) const;
  bool operator != (const Composite& other) const;
};

} // namespace stereopermutation

#endif
