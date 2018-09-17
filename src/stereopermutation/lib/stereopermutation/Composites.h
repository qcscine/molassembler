#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H

#include "chemical_symmetries/Names.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/TinySet.h"
#include "temple/OrderedPair.h"

/*!@file
 *
 * @brief Generate rotational orientations of two symmetries fused over a bond
 *
 * Enables the computation of relative orientations of arbitrary symmetries
 * fused at arbitrary positions that yield indices of permutations that are
 * meaningful for ranking computations.
 */

namespace stereopermutation {

class Composite {
public:
//!@name Member types
//!@{
  struct AngleGroup {
    double angle;
    std::vector<unsigned> symmetryPositions;
    bool isotropic;
  };

  //! Encompasses the orientation of a symmetry along a fused bond
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

    //! Member initializing constructor
    OrientationState(
      Symmetry::Name passSymmetry,
      unsigned passFusedPosition,
      std::vector<char> passCharacters,
      std::size_t passIdentifier
    );

    //! Applies a rotation to the fused position and characters
    void applyCharacterRotation(const std::vector<unsigned>& rotation);

    //! Smallest symmetry position from the same group as the fused position
    unsigned lowestEqualPositionInSymmetry() const;

    //! Calculates the required reduction mapping to the canonical form
    std::vector<unsigned> findReductionMapping(unsigned reducedFusedPosition) const;

    /* c++17 nodiscard */
    /*! Transforms the OrientationState to a canonical form
     *
     * Transforms the OrientationState by applying a reduction mapping to the
     * smallest symmetry position from the same group as the fused position.
     * Returns the mapping needed to revert the OrientationState back to its
     * original data values.
     */
    std::vector<unsigned> transformToCanonical();

    //! Reverts the OrientationState to non-canonical form
    void revert(const std::vector<unsigned>& reversionMapping);

    //! Collects all coplanar indices that are closest to the fused symmetry position
    AngleGroup smallestAngleGroup() const;

    //! Full member lexicographical comparison in order of declaration
    bool operator < (const OrientationState& other) const;
    //! Full member lexicographical comparison in order of declaration
    bool operator == (const OrientationState& other) const;
  };

  //! First symmetry position, second symmetry position, dihedral angle tuple
  using DihedralTuple = std::tuple<unsigned, unsigned, double>;

  //! List of lists of dihedral angles for distinct rotational configurations
  using PermutationsList = std::vector<
    std::vector<DihedralTuple>
  >;

  using PerpendicularAngleGroups = std::vector<
    std::pair<
      temple::TinyUnorderedSet<double>,
      std::vector<
        std::pair<unsigned, unsigned>
      >
    >
  >;
//!@}

public:
//!@name Constructors
//!@{
  Composite(OrientationState first, OrientationState second);
//!@}
//
//!@name Static members
//!@{
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
    Symmetry::Name symmetryName,
    unsigned fixedSymmetryPosition,
    const std::vector<unsigned>& changedPositions
  );

  static std::vector<unsigned> rotation(
    Symmetry::Name symmetryName,
    unsigned fixedSymmetryPosition,
    const std::vector<unsigned>& perpendicularPlanePositions
  );

  //! Creates sets of within-group cross angles in the perpendicular plane
  static PerpendicularAngleGroups inGroupAngles(
    const AngleGroup& angleGroup,
    Symmetry::Name symmetryName
  );
//!@}

//!@name Information
//!@{
  /*! Returns a set of dihedrals for a particular permutation
   *
   * \note The first two elements of each tuple specify the symmetry position
   * within that side's symmetry. The first element is for the left symmetry,
   * the second for the right symmetry.
   */
  const std::vector<DihedralTuple>& dihedrals(unsigned permutationIndex) const;

  //! Returns whether the Composite is isotropic overall
  bool isIsotropic() const;

  //! Returns the orientation state of the composite
  const temple::OrderedPair<OrientationState>& orientations() const;

  //! Returns the number of permutations for this Composite
  unsigned permutations() const;
//!@}

//!@name Iterators
//!@{
  //! Through-iteration of the generated dihedral permutations
  PermutationsList::const_iterator begin() const;
  PermutationsList::const_iterator end() const;
//!@}

//!@name Operators
//!@{
  bool operator == (const Composite& other) const;
  bool operator != (const Composite& other) const;
//!@}

private:
//!@name Private state
//!@{
  //! Stores the relative orientation with which the permutations were generated
  temple::OrderedPair<OrientationState> _orientations;

  //! List of dihedral sets that comprise all spatial arrangements
  PermutationsList _stereopermutations;

  //! Stores whether the Composite is isotropic
  bool _isotropic;
//!@}

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
    double angleFromBoundSymmetryPosition,
    double angleBetweenSubstituents
  );
};

} // namespace stereopermutation

#endif
