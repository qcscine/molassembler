#include "Primitives.h"

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
 */

constexpr decltype(Linear::name) Linear::name;
constexpr decltype(Linear::size) Linear::size;
constexpr decltype(Linear::coordinates) Linear::coordinates;
constexpr decltype(Linear::rotations) Linear::rotations;
constexpr decltype(Linear::tetrahedra) Linear::tetrahedra;
constexpr decltype(Linear::stringName) Linear::stringName;

constexpr decltype(Bent::name) Bent::name;
constexpr decltype(Bent::size) Bent::size;
constexpr decltype(Bent::coordinates) Bent::coordinates;
constexpr decltype(Bent::rotations) Bent::rotations;
constexpr decltype(Bent::tetrahedra) Bent::tetrahedra;
constexpr decltype(Bent::stringName) Bent::stringName;

constexpr decltype(TrigonalPlanar::name) TrigonalPlanar::name;
constexpr decltype(TrigonalPlanar::size) TrigonalPlanar::size;
constexpr decltype(TrigonalPlanar::coordinates) TrigonalPlanar::coordinates;
constexpr decltype(TrigonalPlanar::rotations) TrigonalPlanar::rotations;
constexpr decltype(TrigonalPlanar::tetrahedra) TrigonalPlanar::tetrahedra;
constexpr decltype(TrigonalPlanar::stringName) TrigonalPlanar::stringName;

constexpr decltype(TrigonalPyramidal::name) TrigonalPyramidal::name;
constexpr decltype(TrigonalPyramidal::size) TrigonalPyramidal::size;
constexpr decltype(TrigonalPyramidal::coordinates) TrigonalPyramidal::coordinates;
constexpr decltype(TrigonalPyramidal::rotations) TrigonalPyramidal::rotations;
constexpr decltype(TrigonalPyramidal::tetrahedra) TrigonalPyramidal::tetrahedra;
constexpr decltype(TrigonalPyramidal::stringName) TrigonalPyramidal::stringName;

constexpr decltype(TShaped::name) TShaped::name;
constexpr decltype(TShaped::size) TShaped::size;
constexpr decltype(TShaped::coordinates) TShaped::coordinates;
constexpr decltype(TShaped::rotations) TShaped::rotations;
constexpr decltype(TShaped::tetrahedra) TShaped::tetrahedra;
constexpr decltype(TShaped::stringName) TShaped::stringName;

constexpr decltype(Tetrahedral::name) Tetrahedral::name;
constexpr decltype(Tetrahedral::size) Tetrahedral::size;
constexpr decltype(Tetrahedral::coordinates) Tetrahedral::coordinates;
constexpr decltype(Tetrahedral::rotations) Tetrahedral::rotations;
constexpr decltype(Tetrahedral::tetrahedra) Tetrahedral::tetrahedra;
constexpr decltype(Tetrahedral::stringName) Tetrahedral::stringName;

constexpr decltype(SquarePlanar::name) SquarePlanar::name;
constexpr decltype(SquarePlanar::size) SquarePlanar::size;
constexpr decltype(SquarePlanar::coordinates) SquarePlanar::coordinates;
constexpr decltype(SquarePlanar::rotations) SquarePlanar::rotations;
constexpr decltype(SquarePlanar::tetrahedra) SquarePlanar::tetrahedra;
constexpr decltype(SquarePlanar::stringName) SquarePlanar::stringName;

constexpr decltype(Seesaw::name) Seesaw::name;
constexpr decltype(Seesaw::size) Seesaw::size;
constexpr decltype(Seesaw::coordinates) Seesaw::coordinates;
constexpr decltype(Seesaw::rotations) Seesaw::rotations;
constexpr decltype(Seesaw::tetrahedra) Seesaw::tetrahedra;
constexpr decltype(Seesaw::stringName) Seesaw::stringName;

constexpr decltype(SquarePyramidal::name) SquarePyramidal::name;
constexpr decltype(SquarePyramidal::size) SquarePyramidal::size;
constexpr decltype(SquarePyramidal::coordinates) SquarePyramidal::coordinates;
constexpr decltype(SquarePyramidal::rotations) SquarePyramidal::rotations;
constexpr decltype(SquarePyramidal::tetrahedra) SquarePyramidal::tetrahedra;
constexpr decltype(SquarePyramidal::stringName) SquarePyramidal::stringName;

constexpr decltype(TrigonalBiPyramidal::name) TrigonalBiPyramidal::name;
constexpr decltype(TrigonalBiPyramidal::size) TrigonalBiPyramidal::size;
constexpr decltype(TrigonalBiPyramidal::coordinates) TrigonalBiPyramidal::coordinates;
constexpr decltype(TrigonalBiPyramidal::rotations) TrigonalBiPyramidal::rotations;
constexpr decltype(TrigonalBiPyramidal::tetrahedra) TrigonalBiPyramidal::tetrahedra;
constexpr decltype(TrigonalBiPyramidal::stringName) TrigonalBiPyramidal::stringName;

constexpr decltype(PentagonalPlanar::name) PentagonalPlanar::name;
constexpr decltype(PentagonalPlanar::size) PentagonalPlanar::size;
constexpr decltype(PentagonalPlanar::coordinates) PentagonalPlanar::coordinates;
constexpr decltype(PentagonalPlanar::rotations) PentagonalPlanar::rotations;
constexpr decltype(PentagonalPlanar::tetrahedra) PentagonalPlanar::tetrahedra;
constexpr decltype(PentagonalPlanar::stringName) PentagonalPlanar::stringName;

constexpr decltype(Octahedral::name) Octahedral::name;
constexpr decltype(Octahedral::size) Octahedral::size;
constexpr decltype(Octahedral::coordinates) Octahedral::coordinates;
constexpr decltype(Octahedral::rotations) Octahedral::rotations;
constexpr decltype(Octahedral::tetrahedra) Octahedral::tetrahedra;
constexpr decltype(Octahedral::stringName) Octahedral::stringName;

constexpr decltype(TrigonalPrismatic::name) TrigonalPrismatic::name;
constexpr decltype(TrigonalPrismatic::size) TrigonalPrismatic::size;
constexpr decltype(TrigonalPrismatic::coordinates) TrigonalPrismatic::coordinates;
constexpr decltype(TrigonalPrismatic::rotations) TrigonalPrismatic::rotations;
constexpr decltype(TrigonalPrismatic::tetrahedra) TrigonalPrismatic::tetrahedra;
constexpr decltype(TrigonalPrismatic::stringName) TrigonalPrismatic::stringName;

constexpr decltype(PentagonalPyramidal::name) PentagonalPyramidal::name;
constexpr decltype(PentagonalPyramidal::size) PentagonalPyramidal::size;
constexpr decltype(PentagonalPyramidal::coordinates) PentagonalPyramidal::coordinates;
constexpr decltype(PentagonalPyramidal::rotations) PentagonalPyramidal::rotations;
constexpr decltype(PentagonalPyramidal::tetrahedra) PentagonalPyramidal::tetrahedra;
constexpr decltype(PentagonalPyramidal::stringName) PentagonalPyramidal::stringName;

constexpr decltype(PentagonalBiPyramidal::name) PentagonalBiPyramidal::name;
constexpr decltype(PentagonalBiPyramidal::size) PentagonalBiPyramidal::size;
constexpr decltype(PentagonalBiPyramidal::coordinates) PentagonalBiPyramidal::coordinates;
constexpr decltype(PentagonalBiPyramidal::rotations) PentagonalBiPyramidal::rotations;
constexpr decltype(PentagonalBiPyramidal::tetrahedra) PentagonalBiPyramidal::tetrahedra;
constexpr decltype(PentagonalBiPyramidal::stringName) PentagonalBiPyramidal::stringName;

constexpr decltype(SquareAntiPrismatic::name) SquareAntiPrismatic::name;
constexpr decltype(SquareAntiPrismatic::size) SquareAntiPrismatic::size;
constexpr decltype(SquareAntiPrismatic::coordinates) SquareAntiPrismatic::coordinates;
constexpr decltype(SquareAntiPrismatic::rotations) SquareAntiPrismatic::rotations;
constexpr decltype(SquareAntiPrismatic::tetrahedra) SquareAntiPrismatic::tetrahedra;
constexpr decltype(SquareAntiPrismatic::stringName) SquareAntiPrismatic::stringName;


#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
constexpr decltype(SquareAntiPrismatic::angleLookupTable) SquareAntiPrismatic::angleLookupTable;
#endif

} // namespace data

} // namespace Symmetry
