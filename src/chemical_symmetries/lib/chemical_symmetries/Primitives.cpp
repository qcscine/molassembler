/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "chemical_symmetries/Primitives.h"

#define DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(ShapeClass) \
  constexpr decltype(ShapeClass::shape) ShapeClass::shape; \
  constexpr decltype(ShapeClass::pointGroup) ShapeClass::pointGroup; \
  constexpr decltype(ShapeClass::size) ShapeClass::size; \
  constexpr decltype(ShapeClass::coordinates) ShapeClass::coordinates; \
  constexpr decltype(ShapeClass::rotations) ShapeClass::rotations; \
  constexpr decltype(ShapeClass::tetrahedra) ShapeClass::tetrahedra; \
  constexpr decltype(ShapeClass::stringName) ShapeClass::stringName; \
  constexpr decltype(ShapeClass::mirror) ShapeClass::mirror;

#define DECLARE_CONSTEXPR_ANGLE_LOOKUP(ShapeClass) \
  constexpr decltype(ShapeClass::angleLookupTable) ShapeClass::angleLookupTable;

namespace Scine {

namespace Symmetry {

namespace data {

/* This looks quite heavily as if it still violates DRY, but no idea how to
 * refactor further with the preprocessor
 *
 * NOTE: Although clang-tidy marks the entire block of declarations below as
 * redundant declarations, they are definitely NOT redundant.
 * - Removing the grouped declarations for each Symmetry below, which are all
 *   constexpr, makes it impossible to use these values outside of a constexpr
 *   context.  However, since symmetryData is generated from these classes at
 *   run-time, they need to be available then too.
 * - For the stringName declarations, static const non-literal members require
 *   out-of-class initializers, so they too are non-redundant
 * - Since coordinates, rotations and tetrahedra all have size-specific or
 *   symmetry-specific array lengths, and hence different types, it is better
 *   to just use decltype(...).
 */

DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Line)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Bent)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(EquilateralTriangle)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(VacantTetrahedron)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(T)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Tetrahedron)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Square)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Seesaw)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(TrigonalPyramid)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(SquarePyramid)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(TrigonalBipyramid)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Pentagon)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Octahedron)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(TrigonalPrism)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(PentagonalPyramid)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Hexagon)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(PentagonalBipyramid)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(CappedOctahedron)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(CappedTrigonalPrism)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(SquareAntiprism)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Cube)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(TrigonalDodecahedron)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(HexagonalBipyramid)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(TricappedTrigonalPrism)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(CappedSquareAntiprism)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(HeptagonalBipyramid)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(BicappedSquareAntiprism)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(EdgeContractedIcosahedron)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Icosahedron)
DECLARE_CONSTEXPR_SHAPE_CLASS_MEMBERS(Cuboctahedron)

DECLARE_CONSTEXPR_ANGLE_LOOKUP(TrigonalPrism)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(CappedOctahedron)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(CappedTrigonalPrism)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(SquareAntiprism)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(Cube)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(TrigonalDodecahedron)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(HexagonalBipyramid)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(TricappedTrigonalPrism)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(CappedSquareAntiprism)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(HeptagonalBipyramid)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(BicappedSquareAntiprism)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(EdgeContractedIcosahedron)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(Icosahedron)
DECLARE_CONSTEXPR_ANGLE_LOOKUP(Cuboctahedron)

} // namespace data

} // namespace Symmetry

} // namespace Scine
