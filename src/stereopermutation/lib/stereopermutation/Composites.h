/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Generate rotational orientations of two symmetries fused over a bond
 *
 * Enables the computation of relative orientations of arbitrary symmetries
 * fused at arbitrary positions that yield indices of permutations that are
 * meaningful for ranking computations.
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H

#include "chemical_symmetries/Shapes.h"
#include "temple/OrderedPair.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/OperatorSuppliers.h"

#include <vector>

namespace Scine {

namespace stereopermutation {

/**
 * @brief Represents the composite of two symmetries joined by a bond at
 *   arbitrary symmetry positions
 *
 * Acts as a container-like type for the generated stereopermutations after
 * construction.
 *
 * @todo There's some dead code here I think (see the static functions) and
 * PerpendicularAngleGroups really isn't a telling type
 */
class Composite {
public:
//!@name Member types
//!@{
  //! A group of symmetry positions at an angle from the fused position
  struct AngleGroup {
    //! The angle this group is placed at from the fused position
    double angle;
    //! The symmetry positions making up this group
    std::vector<unsigned> symmetryPositions;
    /*!
     * @brief Whether the ranking characters indicate that this group of
     *   symmetry positions is isotropic
     */
    bool isotropic;
  };

  /*! @brief Encompasses the orientation of a symmetry along a fused bond
   *
   * Comparison is based upon lexicographic comparison of its struct members in
   * order of declaration.
   */
  struct OrientationState : public temple::crtp::LexicographicComparable<OrientationState> {
    //! The symmetry of either positional symmetry
    Symmetry::Shape symmetry;
    //! The symmetry position of the local symmetry that the other is fused at
    unsigned fusedPosition;
    //! Abstract ranking-characters of the ligands at their symmetry positions
    std::vector<char> characters;
    /*!
     * @brief An identifier to the symmetry source
     *
     * Since OrientationState is used internally in an ordered pair which may
     * swap elements on changes, an identifier can be useful in reassociating
     * the OrientationState with its data source.
     */
    std::size_t identifier;

    //! Member initializing constructor
    OrientationState(
      Symmetry::Shape passSymmetry,
      unsigned passFusedPosition,
      std::vector<char> passCharacters,
      std::size_t passIdentifier
    );

    /*! @brief Applies a rotation to the fused position and ranking characters
     *
     * @complexity{@math{\Theta(N)}}
     */
    void applyCharacterRotation(const std::vector<unsigned>& rotation);

    /*! @brief Smallest symmetry position from the same group as the fused position
     *
     * @complexity{@math{\Theta(S^2)}}
     */
    unsigned lowestEqualPositionInSymmetry() const;

    /*! @brief Calculates the required reduction mapping to the canonical form
     *
     * @complexity{@math{O(S!)} where @math{S} is the size of the symmetry.
     * More precisely, this function scales linearly with the number of
     * unlinked stereopermutations if all substituents are different, but that
     * itself probably scales factorially in the size of the symmetry.}
     *
     * @todo Check what this code shares with generateAllRotations and consider
     *   refactoring
     */
    std::vector<unsigned> findReductionMapping(unsigned reducedFusedPosition) const;

    /* c++17 nodiscard */
    /*!
     * @brief Transforms the OrientationState to a canonical form
     *
     * Transforms the OrientationState by applying a reduction mapping to the
     * smallest symmetry position from the same group as the fused position.
     * Returns the mapping needed to revert the OrientationState back to its
     * original data values.
     *
     * @note To transform back to the original state, keep the result of this
     * function and call revert() with it.
     *
     * @complexity{@math{O(S!)} where @math{S} is the size of the symmetry (see
     * findReductionMapping)}
     */
    std::vector<unsigned> transformToCanonical();

    /*! @brief Reverts the OrientationState to a non-canonical form
     *
     * Call this with the result of transformToCanonical()!
     *
     * @complexity{@math{\Theta(S)}}
     */
    void revert(const std::vector<unsigned>& reversionMapping);

    /*!
     * @brief Collects all coplanar indices that are closest to the fused
     *   symmetry position
     *
     * @complexity{@math{\Theta(S)}}
     * @post AngleGroup's symmetry positions are sorted
     */
    AngleGroup smallestAngleGroup() const;

    //! Full tuple-like lexicographical comparison of members in order of declaration
    inline auto tie() const {
      return std::tie(symmetry, fusedPosition, characters);
    }
  };

  //! First symmetry position, second symmetry position, dihedral angle tuple
  using DihedralTuple = std::tuple<unsigned, unsigned, double>;

  //! List of lists of dihedral angles for distinct rotational configurations
  using PermutationsList = std::vector<
    std::vector<DihedralTuple>
  >;

  using PerpendicularAngleGroups = std::vector<
    std::pair<
      std::vector<double>,
      std::vector<
        std::pair<unsigned, unsigned>
      >
    >
  >;

  //! Relative orientation of substituent groups along the dihedral
  enum class Alignment {
    //! At least two substituents eclipse one another along the axis
    Eclipsed,
    //! At least one pair of substituents are staggered along the axis
    Staggered
  };
//!@}

//!@name Constructors
//!@{
  /*! @brief Constructor calculating all permutations
   *
   * @complexity{@math{O(S!)} where S is the size of the larger symmetry of the
   * two OrientationState instances}
   *
   * @post Each permutations' dihedrals are sorted (lexicographically)
   */
  Composite(
    OrientationState first,
    OrientationState second,
    Alignment alignment = Alignment::Eclipsed
  );
//!@}
//
//!@name Static members
//!@{
  //! 1° absolute tolerance floating point comparison helper for angle groups
  static constexpr temple::floating::ExpandedAbsoluteEqualityComparator<double> fpComparator {
    temple::Math::toRadians(1.0)
  };

  /*!
   * @brief Generates a symmetry rotation subject to constraints
   *
   * @param shape The symmetry in which the rotation is sought
   * @param fixedSymmetryPosition The symmetry position in shape that is
   *   to be kept constant
   * @param changedPositions The symmetry positions which must change in the
   *   sought rotation
   */
  static std::vector<unsigned> generateRotation(
    Symmetry::Shape shape,
    unsigned fixedSymmetryPosition,
    const std::vector<unsigned>& changedPositions
  );

  static std::vector<unsigned> rotation(
    Symmetry::Shape shape,
    unsigned fixedSymmetryPosition,
    const std::vector<unsigned>& perpendicularPlanePositions
  );

  //! Creates sets of within-group cross angles in the perpendicular plane
  static PerpendicularAngleGroups inGroupAngles(
    const AngleGroup& angleGroup,
    Symmetry::Shape shape
  );
//!@}

//!@name Modification
//!@{
  void applyIdentifierPermutation(const std::vector<std::size_t>& permutation);
//!@}

//!@name Information
//!@{
  /*!
   * @brief Returns a set of dihedrals for a particular permutation
   *
   * @complexity{@math{\Theta(1)}}
   * @note The first two elements of each tuple specify the symmetry position
   * within that side's symmetry. The first element is for the left symmetry,
   * the second for the right symmetry.
   */
  const std::vector<DihedralTuple>& dihedrals(unsigned permutationIndex) const;

  /*! @brief Returns whether the Composite is isotropic overall
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool isIsotropic() const;

  /*! @brief Returns the higher number of relevant substituent symmetry indices
   *   of both sides
   */
  unsigned order() const;

  /*! @brief Returns the orientation state of the composite
   *
   * @complexity{@math{\Theta(1)}}
   */
  const temple::OrderedPair<OrientationState>& orientations() const;

  /*! @brief Returns the number of permutations for this Composite
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned permutations() const;

  /*! @brief Returns the alignment with which the Composite was generated with
   *
   * @complexity{@math{\Theta(1)}}
   */
  Alignment alignment() const;
//!@}

//!@name Iterators
//!@{
  //! Through-iteration of the generated dihedral permutations
  PermutationsList::const_iterator begin() const;
  PermutationsList::const_iterator end() const;
//!@}

//!@name Operators
//!@{
  //! @brief Less-than comparison by orientations()
  bool operator < (const Composite& other) const;
  //! @brief Equality comparison by orientations()
  bool operator == (const Composite& other) const;
  //! @brief Inequality comparison by orientations()
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

  //! Stores with which Alignment the stereopermutations were generated
  Alignment _alignment;
//!@}

  /*!
   * @brief Calculates the angle between two substituents that have the same
   *   angle from the bound symmetry position in a plane perpendicular to their
   *   shared axis.
   *
   * The axis is defined by the bound symmetry position and the central point of
   * the symmetry. Symmetry positions with identical angles from the bound
   * symmetry position are situated together in a plane perpendicular to the
   * axis. The angle between these symmetry positions via the intersection with
   * the axis in that plane is calculated by this function.
   *
   * @warning Do not call this with angleFromBoundSymmetryPosition == M_PI as
   *   the geometrical idea collapses.
   */
  static double perpendicularSubstituentAngle(
    double angleFromBoundSymmetryPosition,
    double angleBetweenSubstituents
  );
};

} // namespace stereopermutation

} // namespace Scine

#endif
