/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Generate rotational orientations of two symmetries fused over a bond
 *
 * Enables the computation of relative orientations of arbitrary symmetries
 * fused at arbitrary positions that yield indices of permutations that are
 * meaningful for ranking computations.
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H

#include "shapes/Shapes.h"
#include "shapes/Data.h"
#include "temple/OrderedPair.h"
#include "temple/constexpr/FloatingPointComparison.h"

#include <vector>

namespace Scine {
namespace stereopermutation {

/**
 * @brief Represents the composite of two shapes joined by a bond at
 *   arbitrary shape vertices
 *
 * Acts as a container-like type for the generated stereopermutations after
 * construction.
 */
class Composite {
public:
//!@name Member types
//!@{
  //! A group of shape vertices at an angle from the fused position
  struct AngleGroup {
    //! The angle this group is placed at from the fused position
    double angle;
    //! The shape vertices making up this group
    std::vector<shapes::Vertex> vertices;
    /*!
     * @brief Whether the ranking characters indicate that this group of
     *   shape vertices is isotropic
     */
    bool isotropic;
  };

  /*! @brief Encompasses the orientation of a shape along a fused bond
   *
   * Comparison is based upon lexicographic comparison of its struct members in
   * order of declaration.
   */
  struct OrientationState : public temple::crtp::LexicographicComparable<OrientationState> {
    //! The shape of either positional shape
    shapes::Shape shape;
    //! The shape vertex of the local shape that the other is fused at
    shapes::Vertex fusedVertex;
    //! Abstract ranking-characters of the sites at their shape vertices
    std::vector<char> characters;
    /*!
     * @brief An identifier to the shape source
     *
     * Since OrientationState is used internally in an ordered pair which may
     * swap elements on changes, an identifier can be useful in reassociating
     * the OrientationState with its data source.
     */
    std::size_t identifier;

    //! Member initializing constructor
    OrientationState(
      shapes::Shape passShape,
      shapes::Vertex passFusedVertex,
      std::vector<char> passCharacters,
      std::size_t passIdentifier
    );

    /*! @brief Applies a rotation to the fused position and ranking characters
     *
     * @complexity{@math{\Theta(N)}}
     */
    void applyCharacterRotation(const std::vector<shapes::Vertex>& rotation);

    /*! @brief Smallest shape vertex from the same group as the fused position
     *
     * @complexity{@math{\Theta(S^2)}}
     */
    shapes::Vertex lowestEqualVertexInShape() const;

    /*! @brief Calculates the required reduction mapping to the canonical form
     *
     * @complexity{@math{O(S!)} where @math{S} is the size of the shape.
     * More precisely, this function scales linearly with the number of
     * unlinked stereopermutations if all substituents are different, but that
     * itself probably scales factorially in the size of the shape.}
     */
    std::vector<shapes::Vertex> findReductionMapping(shapes::Vertex reducedFusedVertex) const;

    /* c++17 nodiscard */
    /*!
     * @brief Transforms the OrientationState to a canonical form
     *
     * Transforms the OrientationState by applying a reduction mapping to the
     * smallest shape vertex from the same group as the fused position.
     * Returns the mapping needed to revert the OrientationState back to its
     * original data values.
     *
     * @note To transform back to the original state, keep the result of this
     * function and call revert() with it.
     *
     * @complexity{@math{O(S!)} where @math{S} is the size of the shape (see
     * findReductionMapping)}
     */
    std::vector<shapes::Vertex> transformToCanonical();

    /*! @brief Reverts the OrientationState to a non-canonical form
     *
     * Call this with the result of transformToCanonical()!
     *
     * @complexity{@math{\Theta(S)}}
     */
    void revert(const std::vector<shapes::Vertex>& reversionMapping);

    /*!
     * @brief Collects all coplanar indices that are closest to the fused
     *   shape vertex
     *
     * @complexity{@math{\Theta(S)}}
     * @post AngleGroup's shape positions are sorted
     */
    AngleGroup smallestAngleGroup() const;

    //! Full tuple-like lexicographical comparison of members in order of declaration
    inline auto tie() const {
      return std::tie(shape, fusedVertex, characters);
    }
  };

  //! First shape vertex, second shape vertex, dihedral angle tuple
  using DihedralTuple = std::tuple<shapes::Vertex, shapes::Vertex, double>;

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
   * @complexity{@math{O(S!)} where S is the size of the larger shape of the
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
   * @brief Generates a shape rotation subject to constraints
   *
   * @param shape The shape in which the rotation is sought
   * @param fixedVertex The shape vertex in @p shape that is kept constant
   * @param changedVertices The shape vertices which must change in the
   *   sought rotation
   */
  static std::vector<shapes::Vertex> generateRotation(
    shapes::Shape shape,
    shapes::Vertex fixedVertex,
    const std::vector<shapes::Vertex>& changedVertices
  );

  static std::vector<shapes::Vertex> rotation(
    shapes::Shape shape,
    shapes::Vertex fixedVertex,
    const std::vector<shapes::Vertex>& perpendicularPlaneVertices
  );

  //! Creates sets of within-group cross angles in the perpendicular plane
  static PerpendicularAngleGroups inGroupAngles(
    const AngleGroup& angleGroup,
    shapes::Shape shape
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
   * @note The first two elements of each tuple specify the shape position
   * within that side's shape. The first element is for the left shape,
   * the second for the right shape.
   */
  const std::vector<DihedralTuple>& dihedrals(unsigned permutationIndex) const;

  /*! @brief Returns whether the Composite is isotropic overall
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool isIsotropic() const;

  /*! @brief Returns the higher number of relevant substituent shape vertices
   *   of both sides
   *
   * E.g. In a combination of EquilateralTriangle and Tetrahedron, the order is
   * four.
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
  temple::OrderedPair<OrientationState> orientations_;

  //! List of dihedral sets that comprise all spatial arrangements
  PermutationsList stereopermutations_;

  //! Stores whether the Composite is isotropic
  bool isotropic_;

  //! Stores with which Alignment the stereopermutations were generated
  Alignment alignment_;
//!@}

  /*!
   * @brief Calculates the angle between two substituents that have the same
   *   angle from the bound shape vertices in a plane perpendicular to their
   *   shared axis.
   *
   * The axis is defined by the bound shape position and the central point of
   * the shape. shape positions with identical angles from the bound
   * shape position are situated together in a plane perpendicular to the
   * axis. The angle between these shape positions via the intersection with
   * the axis in that plane is calculated by this function.
   *
   * @warning Do not call this with angleFromBoundShapeVertex == M_PI as
   *   the geometrical idea collapses.
   */
  static double perpendicularSubstituentAngle(
    double angleFromBoundShapeVertex,
    double angleBetweenSubstituents
  );
};

} // namespace stereopermutation
} // namespace Scine

#endif
