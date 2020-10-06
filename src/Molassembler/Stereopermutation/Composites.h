/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Generate rotational orientations of two shapes fused over a bond
 *
 * Enables the computation of relative orientations of arbitrary shapes
 * fused at arbitrary positions that yield indices of permutations that are
 * meaningful for ranking computations.
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATION_COMPOSITES_H

#include "Molassembler/Shapes/Shapes.h"
#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Temple/OrderedPair.h"
#include "Molassembler/Temple/constexpr/FloatingPointComparison.h"

#include "boost/optional.hpp"

namespace Scine {
namespace Molassembler {
namespace Stereopermutations {

/**
 * @brief Represents the composite of two shapes joined by a bond at
 *   arbitrary shape vertices
 *
 * Acts as a container-like type for the generated stereopermutations after
 * construction.
 */
class MASM_EXPORT Composite {
public:
//!@name Member types
//!@{
  //! A group of shape vertices at an angle from the fused position
  struct AngleGroup {
    //! The angle this group is placed at from the fused position
    double angle;
    //! The shape vertices making up this group
    std::vector<Shapes::Vertex> vertices;
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
  struct OrientationState : public Temple::Crtp::LexicographicComparable<OrientationState> {
    //! The shape of either positional shape
    Shapes::Shape shape;
    //! The shape vertex of the local shape that the other is fused at
    Shapes::Vertex fusedVertex;
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
      Shapes::Shape passShape,
      Shapes::Vertex passFusedVertex,
      std::vector<char> passCharacters,
      std::size_t passIdentifier
    );

    /*! @brief Applies a rotation to the fused position and ranking characters
     *
     * @complexity{@math{\Theta(N)}}
     */
    std::vector<char> applyCharacterRotation(const std::vector<Shapes::Vertex>& rotation) const;

    /*! @brief Smallest shape vertex from the same group as the fused position
     *
     * @complexity{@math{\Theta(S^2)}}
     */
    Shapes::Vertex lowestEqualVertexInShape() const;

    /*! @brief Calculates the required reduction mapping to the canonical form
     *
     * @complexity{@math{O(S!)} where @math{S} is the size of the shape.
     * More precisely, this function scales linearly with the number of
     * unlinked stereopermutations if all substituents are different, but that
     * itself probably scales factorially in the size of the shape.}
     */
    std::vector<Shapes::Vertex> findReductionMapping(Shapes::Vertex reducedFusedVertex) const;

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
    std::vector<Shapes::Vertex> transformToCanonical();

    /*! @brief Reverts the OrientationState to a non-canonical form
     *
     * Call this with the result of transformToCanonical()!
     *
     * @complexity{@math{\Theta(S)}}
     */
    void revert(const std::vector<Shapes::Vertex>& reversionMapping);

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
      // Comparisons should be performed in canonical space
      const Shapes::Vertex reducedVertex = lowestEqualVertexInShape();
      return std::make_tuple(
        shape,
        reducedVertex,
        applyCharacterRotation(findReductionMapping(reducedVertex))
      );
    }
  };

  //! Relative orientation of substituent groups along the dihedral
  enum class Alignment {
    //! At least two substituents eclipse one another along the axis
    Eclipsed,
    //! At least one pair of substituents are staggered along the axis
    Staggered,
    //! Both eclipsed and staggered alignments
    EclipsedAndStaggered,
    //! Offset exactly halfway between eclipsed and staggered alignments
    BetweenEclipsedAndStaggered
  };

  //! Class representing a single rotational permutation
  struct Permutation {
    using DihedralTuple = std::tuple<Shapes::Vertex, Shapes::Vertex, double>;
    using VertexPair = std::pair<Shapes::Vertex, Shapes::Vertex>;

    bool close(const std::vector<DihedralTuple>& comparisonDihedrals) const;

    //! Shape vertices aligned for this set of dihedrals
    VertexPair alignedVertices;
    //! Set of dihedrals between all relevant vertices along the bond
    std::vector<DihedralTuple> dihedrals;
    //! Shape vertices this permutation is ranking equivalent to, if applicable
    boost::optional<VertexPair> rankingEquivalentTo;
  };

  using PermutationsList = std::vector<Permutation>;

  struct PermutationGenerator {
    static bool isDuplicate(
      Permutation permutation,
      const PermutationsList& permutations,
      double degrees
    );

    static PermutationsList deduplicate(PermutationsList permutations, double degrees);

    PermutationGenerator(Temple::OrderedPair<OrientationState> orientations);

    double dihedral(
      const Shapes::Vertex firstVertex,
      const Shapes::Vertex secondVertex
    ) const;

    Permutation align(
      const Shapes::Vertex firstVertex,
      const Shapes::Vertex secondVertex,
      Alignment alignment
    );

    PermutationsList generateEclipsedOrStaggered(
      Alignment alignment,
      double deduplicationDegrees
    );

    PermutationsList generate(
      Alignment alignment,
      double deduplicationDegrees=15
    );

  //!@name Members
  //!@{
    Temple::OrderedPair<OrientationState> orientations;
    std::pair<Shapes::Permutation, Shapes::Permutation> reversionMappings;
    std::pair<AngleGroup, AngleGroup> angleGroups;
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> coordinates;
  //!@}
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
  static constexpr Temple::Floating::ExpandedAbsoluteEqualityComparator<double> fpComparator {
    Temple::Math::toRadians(1.0)
  };

  /*!
   * @brief Generates a shape rotation subject to constraints
   *
   * @param shape The shape in which the rotation is sought
   * @param fixedVertex The shape vertex in @p shape that is kept constant
   * @param changedVertices The shape vertices which must change in the
   *   sought rotation
   */
  static std::vector<Shapes::Vertex> generateRotation(
    Shapes::Shape shape,
    Shapes::Vertex fixedVertex,
    const std::vector<Shapes::Vertex>& changedVertices
  );

  static std::vector<Shapes::Vertex> rotation(
    Shapes::Shape shape,
    Shapes::Vertex fixedVertex,
    const std::vector<Shapes::Vertex>& perpendicularPlaneVertices
  );
//!@}

//!@name Modification
//!@{
  void applyIdentifierPermutation(const std::vector<std::size_t>& permutation);
//!@}

//!@name Information
//!@{
  //! Returns all permutations, including ranking-duplicates
  const PermutationsList& allPermutations() const;

  /*! @brief Returns whether the Composite is isotropic overall
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool isIsotropic() const;

  /*! @brief Returns the higher number of relevant substituent shape vertices
   *   of both sides
   *
   * E.g. In a combination of EquilateralTriangle and Tetrahedron, the order is
   * three.
   */
  unsigned order() const;

  /*! @brief Returns the orientation state of the composite
   *
   * @complexity{@math{\Theta(1)}}
   */
  const Temple::OrderedPair<OrientationState>& orientations() const;

  //! Index of the ranking equivalent base to a permutation
  unsigned rankingEquivalentBase(const unsigned permutation) const;

  //! List of all base permutation indices (without ranking-equivalents)
  std::vector<unsigned> nonEquivalentPermutationIndices() const;

  /*! @brief Returns the number of non-ranking-equivalent permutations for this
   *   Composite
   *
   * @complexity{@math{\Theta(1)}}
   */
  unsigned countNonEquivalentPermutations() const;

  /*! @brief Returns the alignment with which the Composite was generated with
   *
   * @complexity{@math{\Theta(1)}}
   */
  Alignment alignment() const;
//!@}

//!@name Iterators
//!@{
  //! Through-iteration of all generated dihedral permutations
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
  Temple::OrderedPair<OrientationState> orientations_;

  //! List of dihedral sets that comprise all spatial arrangements
  PermutationsList stereopermutations_;

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

} // namespace Stereopermutations
} // namespace Molassembler
} // namespace Scine

#endif
