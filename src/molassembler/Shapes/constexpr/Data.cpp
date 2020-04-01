/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Shapes/constexpr/Data.h"

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
namespace shapes {
namespace concepts {

/**
 * @brief Concept checking class for symmetry classes
 *
 * @tparam T Class to check the concept for
 *
 * All members must be `static constexpr` (`const` is then implicit)
 *
 * @par Member requirements (Type and name)
 * - `Shape shape`: Enum member specifying the symmetry name
 * - `unsigned size`: Number of symmetry positions in the symmetry
 * - `char* stringName`: Human readable string of the name
 * - `double(const unsigned, const unsigned) angleFunction`: Angle in radians
 *   between symmetry position indices
 * - `std::array<temple::Vector, size> coordinate`: Origin-centered normalized
 *   position vectors of the symmetry positions
 * - `std::array< std::array<unsigned, size>, ?> rotations`: Spatial rotations
 *   represented as index permutations between symmetry positions. A minimal
 *   set that combined can generate all rotations.
 * - `std::array< std::array<unsigned, 4>, ?> tetrahedra`: A list of tetrahedra
 *   definitions whose signed volumes can uniquely identify a maximally
 *   asymmetric set of ligands
 * - `std::array<unsigned, 0 or size> mirror`: A mirroring symmetry element or
 *   an empty array if the symmetry cannot be chiral.
 */
template<typename T>
struct ShapeClass : std::integral_constant<bool,
  (
    std::is_same<Shape, std::decay_t<decltype(T::shape)>>::value
    && std::is_same<PointGroup, std::decay_t<decltype(T::pointGroup)>>::value
    && std::is_same<unsigned, std::decay_t<decltype(T::size)>>::value
    && std::is_same<const char*, std::decay_t<decltype(T::stringName)>>::value
    && std::is_same<double, decltype(T::angleFunction(0u, 1u))>::value
    && std::is_same<
      temple::Vector,
      temple::getValueType<decltype(T::coordinates)>
    >::value
    && T::coordinates.size() == T::size
    && std::is_same<
      std::array<unsigned, T::size>,
      temple::getValueType<decltype(T::rotations)>
    >::value
    && std::is_same<
      std::array<unsigned, 4>,
      temple::getValueType<decltype(T::tetrahedra)>
    >::value
    && std::is_same<unsigned, temple::getValueType<decltype(T::mirror)>>::value
    && (T::mirror.size() == 0 || T::mirror.size() == T::size)
  )
> {};

template<typename T>
constexpr bool isIotaPermutation(const T& indices) {
  const unsigned size = indices.size();

  for(unsigned i = 0; i < size; ++i) {
    bool foundI = false;

    for(unsigned j = 0; j < size; ++j) {
      if(indices[j] == i) {
        foundI = true;
        break;
      }
    }

    if(!foundI) {
      return false;
    }
  }

  return true;
}

template<typename T>
constexpr bool all_of(const T& t) {
  for(unsigned i = 0; i < t.size(); ++i) {
    if(!t[i]) {
      return false;
    }
  }

  return true;
}

template<typename T, std::size_t ... Inds>
constexpr bool allRotationsValid(std::index_sequence<Inds ...> /* inds */) {
  constexpr unsigned numRotations = sizeof...(Inds);
  temple::Array<bool, numRotations> valid {
    isIotaPermutation(T::rotations[Inds])...
  };
  return all_of(valid);
}

template<typename T>
struct ValidRotations : std::integral_constant<bool,
  allRotationsValid<T>(
    std::make_index_sequence<T::rotations.size()> {}
  )
> {};

template<typename T>
struct ValidMirror : std::integral_constant<bool,
  isIotaPermutation(T::mirror)
> {};

template<typename T, std::size_t ... Inds>
constexpr bool allVectorsNormalized(std::index_sequence<Inds ...> /* inds */) {
  temple::Array<bool, T::size> valid {
    (temple::Math::abs(T::coordinates[Inds].norm() - 1) < 1e-6)...
  };
  return all_of(valid);
}

template<typename T>
struct ValidCoordinates : std::integral_constant<bool,
  allVectorsNormalized<T>(std::make_index_sequence<T::size> {})
> {};

} // namespace concepts

namespace data {

/* Static property correctness checking */

static_assert(
  temple::tuples::allOf<allShapeDataTypes, concepts::ShapeClass>(),
  "Not all shape data types fulfill the ShapeClass concept"
);

static_assert(
  std::tuple_size<allShapeDataTypes>::value == nShapes,
  "Not all shape names have a shape type"
);

static_assert(
  temple::tuples::allOf<allShapeDataTypes, concepts::ValidRotations>(),
  "Not all shape data types' rotations are valid"
);

static_assert(
  temple::tuples::allOf<allShapeDataTypes, concepts::ValidMirror>(),
  "Not all shape data types' mirrors are valid"
);

static_assert(
  temple::tuples::allOf<allShapeDataTypes, concepts::ValidCoordinates>(),
  "Not all shape data types' coordinates are valid"
);

/* This looks quite heavily as if it still violates DRY, but no idea how to
 * refactor further with the preprocessor
 *
 * NOTE: Although clang-tidy marks the entire block of declarations below as
 * redundant declarations, they are definitely NOT redundant.
 * - Removing the grouped declarations for each Symmetry below, which are all
 *   constexpr, makes it impossible to use these values outside of a constexpr
 *   context.  However, since shapeData is generated from these classes at
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
} // namespace shapes
} // namespace Scine
