/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Interface for shape property caching
 */

#ifndef INCLUDE_SHAPE_PROPERTY_CACHING_H
#define INCLUDE_SHAPE_PROPERTY_CACHING_H

#include "molassembler/Shapes/constexpr/Properties.h"
#include "molassembler/Shapes/Properties.h"

namespace Scine {
namespace Shapes {

//! Precomputed min and max angle values in radians for all symmetries
extern const Temple::Array<std::pair<double, double>, nShapes> symmetryAngleBounds;

/*! @brief Calculate the minimum angle in a symmetry
 *
 * Calculates the minimum angle between symmetry positions in a symmetry class
 *
 * @complexity{@math{\Theta(S^2)}}
 * @see ConstexprProperties::calculateSmallestAngle
 */
MASM_EXPORT double minimumAngle(Shape symmetryName);

/*! @brief Calculate the maximum angle in a symmetry
 *
 * Calculates the maximum angle between symmetry positions in a symmetry class
 *
 * @complexity{@math{\Theta(S^2)}}
 */
MASM_EXPORT double maximumAngle(Shape symmetryName);

/* Derived stored constexpr data */
/*! @brief The smallest angle between ligands in all symmetries
 *
 * Stores the smallest angle between symmetry positions across all symmetries
 *
 * @complexity{@math{\Theta(NS^2)} where @math{N} is the number of largest symmetries and @math{S} is the size of the largest symmetry}
 */
constexpr double smallestAngle [[gnu::unused]]
= Temple::Tuples::unpackToFunction<
  Data::allShapeDataTypes,
  ConstexprProperties::minAngleFunctor
>();

#ifdef USE_CONSTEXPR_TRANSITION_MAPPINGS
/*! @brief 0, +1 symmetry transition mappings calculated at compile-time
 *
 * A strictly upper triangular matrix of mapping optionals. Only transitions
 * between symmetries of equivalent or increasing sizes are populated, all
 * other are Nones.
 *
 * @complexity{@math{\Theta(S!)} where @math{S} is the size of the largest symmetry}
 */
extern const Temple::UpperTriangularMatrix<
  Temple::Optional<ConstexprProperties::MappingsReturnType>,
  nShapes * (nShapes - 1) / 2
> allMappings;
#endif

/* Dynamic access to constexpr data */
/*! @brief Cache for on-the-fly generated mappings between symmetries
 *
 * Accesses allMappings if it was generated.
 */
extern Temple::MinimalCache<
  std::tuple<Shape, Shape, boost::optional<unsigned>>,
  Properties::ShapeTransitionGroup
> mappingsCache;

/*! @brief Cached access to mappings. Populates the cache from constexpr if generated.
 *
 * @complexity{@math{\Theta(S!)} where @math{S} is the size of the symmetry if
 * the transition is not cached, @math{\Theta(1)} otherwise.}
 *
 * @param a Source shape
 * @param b Target shape
 * @param removedIndexOption The symmetry position removed from @p source if the
 *   transition is to be a symmetry position loss. Necessary if the size of the
 *   target is one less than that of the source. Defaults to None.
 *
 * @returns The symmetry transition if possible, None otherwise
 */
const boost::optional<const Properties::ShapeTransitionGroup&> getMapping(
  Shape a,
  Shape b,
  const boost::optional<Vertex>& removedIndexOption = boost::none
);

#ifdef USE_CONSTEXPR_HAS_MULTIPLE_UNLINKED_STEREOPERMUTATIONS
/*! @brief All precomputed values for hasMultipleUnlinkedStereopermutations
 *
 * Precomputes the value of hasMultipleUnlinkedStereopermutations for all
 * symmetries if enabled.
 *
 * @complexity{@math{\Theta(S!)} where @math{S} is the size of the largest symmetry}
 */
extern const Temple::Array<
  Temple::DynamicArray<bool, ConstexprProperties::maxShapeSize>,
  nShapes
> allHasMultipleUnlinkedStereopermutations;
#endif

//! Run-time cache
extern Temple::MinimalCache<
  Shape,
  std::vector<bool>
> hasMultipleUnlinkedCache;

/*! @brief Cached access to multiple unlinked values
 *
 * Populates the cache with allHasMultipleUnlinkedStereopermutations if the
 * constexpr number of unlinked stereopermutations was calculated at runtime.
 * If not, calculates the value and stores it in the cache.
 *
 * @complexity{@math{\Theta(S!)} where @math{S} is the size of the symmetry if
 * the symmetry and number of identical ligands is not cached, @math{\Theta(1)} otherwise}
 *
 * @param symmetryName symmetry for which to check the value
 * @param nIdenticalLigands number of equivalent ligands
 *
 * @returns Whether there are multiple stereopermutations assuming no ligands
 *   are linked
 */
MASM_EXPORT bool hasMultipleUnlinkedStereopermutations(
  Shape symmetryName,
  unsigned nIdenticalLigands
);

} // namespace Shapes
} // namespace Scine

#endif
