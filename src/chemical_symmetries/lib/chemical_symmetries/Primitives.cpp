/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "chemical_symmetries/Primitives.h"

namespace Scine {

namespace Symmetry {

namespace data {

/* This looks quite heavily as if it violates DRY, but no idea how to refactor,
 * not even with preprocessor
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

constexpr decltype(Line::shape) Line::shape;
constexpr decltype(Line::size) Line::size;
constexpr decltype(Line::coordinates) Line::coordinates;
constexpr decltype(Line::rotations) Line::rotations;
constexpr decltype(Line::tetrahedra) Line::tetrahedra;
constexpr decltype(Line::stringName) Line::stringName;
constexpr decltype(Line::mirror) Line::mirror;

constexpr decltype(Bent::shape) Bent::shape;
constexpr decltype(Bent::size) Bent::size;
constexpr decltype(Bent::coordinates) Bent::coordinates;
constexpr decltype(Bent::rotations) Bent::rotations;
constexpr decltype(Bent::tetrahedra) Bent::tetrahedra;
constexpr decltype(Bent::stringName) Bent::stringName;
constexpr decltype(Bent::mirror) Bent::mirror;

constexpr decltype(EquilateralTriangle::shape) EquilateralTriangle::shape;
constexpr decltype(EquilateralTriangle::size) EquilateralTriangle::size;
constexpr decltype(EquilateralTriangle::coordinates) EquilateralTriangle::coordinates;
constexpr decltype(EquilateralTriangle::rotations) EquilateralTriangle::rotations;
constexpr decltype(EquilateralTriangle::tetrahedra) EquilateralTriangle::tetrahedra;
constexpr decltype(EquilateralTriangle::stringName) EquilateralTriangle::stringName;
constexpr decltype(EquilateralTriangle::mirror) EquilateralTriangle::mirror;

constexpr decltype(VacantTetrahedron::shape) ApicalTrigonalPyramid::shape;
constexpr decltype(VacantTetrahedron::size) ApicalTrigonalPyramid::size;
constexpr decltype(VacantTetrahedron::coordinates) ApicalTrigonalPyramid::coordinates;
constexpr decltype(VacantTetrahedron::rotations) ApicalTrigonalPyramid::rotations;
constexpr decltype(VacantTetrahedron::tetrahedra) ApicalTrigonalPyramid::tetrahedra;
constexpr decltype(VacantTetrahedron::stringName) ApicalTrigonalPyramid::stringName;
constexpr decltype(VacantTetrahedron::mirror) ApicalTrigonalPyramid::mirror;

constexpr decltype(T::shape) T::shape;
constexpr decltype(T::size) T::size;
constexpr decltype(T::coordinates) T::coordinates;
constexpr decltype(T::rotations) T::rotations;
constexpr decltype(T::tetrahedra) T::tetrahedra;
constexpr decltype(T::stringName) T::stringName;
constexpr decltype(T::mirror) T::mirror;

constexpr decltype(Tetrahedron::shape) Tetrahedron::shape;
constexpr decltype(Tetrahedron::size) Tetrahedron::size;
constexpr decltype(Tetrahedron::coordinates) Tetrahedron::coordinates;
constexpr decltype(Tetrahedron::rotations) Tetrahedron::rotations;
constexpr decltype(Tetrahedron::tetrahedra) Tetrahedron::tetrahedra;
constexpr decltype(Tetrahedron::stringName) Tetrahedron::stringName;
constexpr decltype(Tetrahedron::mirror) Tetrahedron::mirror;

constexpr decltype(Square::shape) Square::shape;
constexpr decltype(Square::size) Square::size;
constexpr decltype(Square::coordinates) Square::coordinates;
constexpr decltype(Square::rotations) Square::rotations;
constexpr decltype(Square::tetrahedra) Square::tetrahedra;
constexpr decltype(Square::stringName) Square::stringName;
constexpr decltype(Square::mirror) Square::mirror;

constexpr decltype(Seesaw::shape) Disphenoid::shape;
constexpr decltype(Seesaw::size) Disphenoid::size;
constexpr decltype(Seesaw::coordinates) Disphenoid::coordinates;
constexpr decltype(Seesaw::rotations) Disphenoid::rotations;
constexpr decltype(Seesaw::tetrahedra) Disphenoid::tetrahedra;
constexpr decltype(Seesaw::stringName) Disphenoid::stringName;
constexpr decltype(Seesaw::mirror) Disphenoid::mirror;

constexpr decltype(TrigonalPyramid::shape) TrigonalPyramid::shape;
constexpr decltype(TrigonalPyramid::size) TrigonalPyramid::size;
constexpr decltype(TrigonalPyramid::coordinates) TrigonalPyramid::coordinates;
constexpr decltype(TrigonalPyramid::rotations) TrigonalPyramid::rotations;
constexpr decltype(TrigonalPyramid::tetrahedra) TrigonalPyramid::tetrahedra;
constexpr decltype(TrigonalPyramid::stringName) TrigonalPyramid::stringName;
constexpr decltype(TrigonalPyramid::mirror) TrigonalPyramid::mirror;

constexpr decltype(SquarePyramid::shape) SquarePyramid::shape;
constexpr decltype(SquarePyramid::size) SquarePyramid::size;
constexpr decltype(SquarePyramid::coordinates) SquarePyramid::coordinates;
constexpr decltype(SquarePyramid::rotations) SquarePyramid::rotations;
constexpr decltype(SquarePyramid::tetrahedra) SquarePyramid::tetrahedra;
constexpr decltype(SquarePyramid::stringName) SquarePyramid::stringName;
constexpr decltype(SquarePyramid::mirror) SquarePyramid::mirror;

constexpr decltype(TrigonalBipyramid::shape) TrigonalBipyramid::shape;
constexpr decltype(TrigonalBipyramid::size) TrigonalBipyramid::size;
constexpr decltype(TrigonalBipyramid::coordinates) TrigonalBipyramid::coordinates;
constexpr decltype(TrigonalBipyramid::rotations) TrigonalBipyramid::rotations;
constexpr decltype(TrigonalBipyramid::tetrahedra) TrigonalBipyramid::tetrahedra;
constexpr decltype(TrigonalBipyramid::stringName) TrigonalBipyramid::stringName;
constexpr decltype(TrigonalBipyramid::mirror) TrigonalBipyramid::mirror;

constexpr decltype(Pentagon::shape) Pentagon::shape;
constexpr decltype(Pentagon::size) Pentagon::size;
constexpr decltype(Pentagon::coordinates) Pentagon::coordinates;
constexpr decltype(Pentagon::rotations) Pentagon::rotations;
constexpr decltype(Pentagon::tetrahedra) Pentagon::tetrahedra;
constexpr decltype(Pentagon::stringName) Pentagon::stringName;
constexpr decltype(Pentagon::mirror) Pentagon::mirror;

constexpr decltype(Octahedron::shape) Octahedron::shape;
constexpr decltype(Octahedron::size) Octahedron::size;
constexpr decltype(Octahedron::coordinates) Octahedron::coordinates;
constexpr decltype(Octahedron::rotations) Octahedron::rotations;
constexpr decltype(Octahedron::tetrahedra) Octahedron::tetrahedra;
constexpr decltype(Octahedron::stringName) Octahedron::stringName;
constexpr decltype(Octahedron::mirror) Octahedron::mirror;

constexpr decltype(TrigonalPrism::shape) TrigonalPrism::shape;
constexpr decltype(TrigonalPrism::size) TrigonalPrism::size;
constexpr decltype(TrigonalPrism::coordinates) TrigonalPrism::coordinates;
constexpr decltype(TrigonalPrism::rotations) TrigonalPrism::rotations;
constexpr decltype(TrigonalPrism::tetrahedra) TrigonalPrism::tetrahedra;
constexpr decltype(TrigonalPrism::stringName) TrigonalPrism::stringName;
constexpr decltype(TrigonalPrism::mirror) TrigonalPrism::mirror;

constexpr decltype(PentagonalPyramid::shape) PentagonalPyramid::shape;
constexpr decltype(PentagonalPyramid::size) PentagonalPyramid::size;
constexpr decltype(PentagonalPyramid::coordinates) PentagonalPyramid::coordinates;
constexpr decltype(PentagonalPyramid::rotations) PentagonalPyramid::rotations;
constexpr decltype(PentagonalPyramid::tetrahedra) PentagonalPyramid::tetrahedra;
constexpr decltype(PentagonalPyramid::stringName) PentagonalPyramid::stringName;
constexpr decltype(PentagonalPyramid::mirror) PentagonalPyramid::mirror;

constexpr decltype(PentagonalBipyramid::shape) PentagonalBipyramid::shape;
constexpr decltype(PentagonalBipyramid::size) PentagonalBipyramid::size;
constexpr decltype(PentagonalBipyramid::coordinates) PentagonalBipyramid::coordinates;
constexpr decltype(PentagonalBipyramid::rotations) PentagonalBipyramid::rotations;
constexpr decltype(PentagonalBipyramid::tetrahedra) PentagonalBipyramid::tetrahedra;
constexpr decltype(PentagonalBipyramid::stringName) PentagonalBipyramid::stringName;
constexpr decltype(PentagonalBipyramid::mirror) PentagonalBipyramid::mirror;

constexpr decltype(SquareAntiprism::shape) SquareAntiprism::shape;
constexpr decltype(SquareAntiprism::size) SquareAntiprism::size;
constexpr decltype(SquareAntiprism::coordinates) SquareAntiprism::coordinates;
constexpr decltype(SquareAntiprism::rotations) SquareAntiprism::rotations;
constexpr decltype(SquareAntiprism::tetrahedra) SquareAntiprism::tetrahedra;
constexpr decltype(SquareAntiprism::stringName) SquareAntiprism::stringName;
constexpr decltype(SquareAntiprism::mirror) SquareAntiprism::mirror;

#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
constexpr decltype(SquareAntiprism::angleLookupTable) SquareAntiprism::angleLookupTable;
#endif

} // namespace data

} // namespace Symmetry

} // namespace Scine
