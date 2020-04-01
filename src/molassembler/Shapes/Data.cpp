/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "shapes/Data.h"

#include "shapes/constexpr/Data.h"
#include "boost/optional.hpp"
#include "temple/constexpr/TupleType.h"

#include <unordered_map>

namespace Scine {
namespace shapes {

//! Dynamic shape information data struct
struct ShapeInformation {
  const std::string stringName;
  const unsigned size;
  const RotationsList rotations;
  const TetrahedronList tetrahedra;
  const Coordinates coordinates;
  const Permutation mirror;
  const PointGroup pointGroup;
  bool threeDimensional;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

//! Map type used to store shape information structs
using ShapeDataMapType = std::unordered_map<
  Shape,
  ShapeInformation,
  std::hash<Shape>,
  std::equal_to<Shape>,
  Eigen::aligned_allocator<std::pair<const Shape, ShapeInformation>>
>;

namespace data {

Eigen::Vector3d toEigen(const temple::Vector& cVector) {
  return {
    cVector.data[0],
    cVector.data[1],
    cVector.data[2]
  };
}

// Typedef to avoid reusing C-Style function ptr type
using AngleFunctionPtr = double(*)(const unsigned, const unsigned);

/*!
 * Constructs an array of function pointers to all static angle functions
 * for runtime lookup
 */
template<typename ... ShapeClasses>
struct angleFunctionFunctor {
  static constexpr std::array<
    data::AngleFunctionPtr,
    sizeof...(ShapeClasses)
  > value() {
    return {{
      &ShapeClasses::angleFunction...
    }};
  }
};

/*! Conversion function to make the dynamic rotations list type from the
 * constexpr data types given in a specifc shape class type
 */
template<size_t shapeSize, size_t nRotations>
RotationsList makeRotations(
  const std::array<
    std::array<unsigned, shapeSize>,
    nRotations
  >& constexprRotations
) {
  RotationsList rotations;

  for(const auto& rotation : constexprRotations) {
    rotations.emplace_back(
      rotation.begin(),
      rotation.end()
    );
  }

  return rotations;
}

/*! Conversion function to make the dynamic tetrahedron list type from the
 * constexpr data types given in a specifc shape class type
 */
template<size_t nTetrahedra>
TetrahedronList makeTetrahedra(
  const std::array<
    std::array<unsigned, 4>,
    nTetrahedra
  >& constexprTetrahedra
) {
  TetrahedronList tetrahedra (nTetrahedra);

  for(unsigned i = 0; i < nTetrahedra; ++i) {
    for(unsigned j = 0; j < 4; ++j) {
      const auto& constexprValue = constexprTetrahedra.at(i).at(j);

      if(constexprValue == ORIGIN_PLACEHOLDER) {
        tetrahedra.at(i).at(j) = boost::none;
      } else {
        tetrahedra.at(i).at(j) = Vertex(constexprValue);
      }
    }
  }

  return tetrahedra;
}

//! Conversion helper to Eigen type from constexpr vector type
Eigen::Vector3d toEigen(const temple::Vector& cVector);

/*! Conversion function to make the dynamic coordinates list type from the
 * constexpr data types given in a specific shape class type
 */
template<size_t shapeSize>
Coordinates makeCoordinates(
  const std::array<temple::Vector, shapeSize>& constexprCoordinates
) {
  Coordinates coordinates(3, shapeSize);

  for(unsigned i = 0; i < shapeSize; ++i) {
    coordinates.col(i) = toEigen(constexprCoordinates.at(i));
  }

  return coordinates;
}

template<size_t shapeSize>
Permutation makeMirror(
  const std::array<unsigned, shapeSize>& constexprMirror
) {
  std::vector<Vertex> mirror (shapeSize);
  std::copy(
    std::begin(constexprMirror),
    std::end(constexprMirror),
    std::begin(mirror)
  );
  return mirror;
}

/* Figure out whether a shape is three dimensional or not
 *
 * Linear in the size of the shape
 */
template<typename ShapeClass>
constexpr bool isThreeDimensional() {
  if(ShapeClass::size == 2) {
    return false;
  }

  // temple::Vector a {}; // zero-vector
  const temple::Vector& c = ShapeClass::coordinates.at(0);
  const temple::Vector& d = ShapeClass::coordinates.at(1);
  const temple::Vector cMinusD = c - d;

  // All points are within a plane containing the origin and the first two vertices
  for(unsigned i = 2; i < ShapeClass::size; ++i) {
    const temple::Vector& b = ShapeClass::coordinates.at(i);
    /* The first line of this would read a-d, but a is always zero
     *
     * And since we're taking the absolute value of the whole thing anyway, we
     * can reduce a - d to just d, inverting the volume, but not changing the
     * absolute value
     */
    const double volume = d.dot((b - d).cross(cMinusD));

    if(std::fabs(volume) > 1e-5) {
      return true;
    }
  }

  return false;
}

/*! @brief Constructs ShapeInformation instance for a shape class
 *
 * @tparam ShapeClass model of concepts::ShapeClass
 */
template<typename ShapeClass>
ShapeInformation makeShapeInformation() {
  return ShapeInformation {
    ShapeClass::stringName,
    ShapeClass::size,
    makeRotations(ShapeClass::rotations),
    makeTetrahedra(ShapeClass::tetrahedra),
    makeCoordinates(ShapeClass::coordinates),
    makeMirror(ShapeClass::mirror),
    ShapeClass::pointGroup,
    isThreeDimensional<ShapeClass>()
  };
}

/*! @brief Creates a map initialization pair for a specific shape class
 *
 * The key is the name, the mapped_type a ShapeInformation instance
 *
 * @tparam ShapeClass model of concepts::ShapeClass
 */
template<typename ShapeClass>
std::pair<Shape, ShapeInformation> makeMapInitPair() {
  return {
    ShapeClass::shape,
    makeShapeInformation<ShapeClass>()
  };
}

/*! Creates the mapping between a shape class's name and its dynamic
 * information in order to have runtime lookup based on shape names.
 */
template<typename ... ShapeClasses>
struct shapeInformationFunctor {
  static ShapeDataMapType value() {
    return {{
      makeMapInitPair<ShapeClasses>()...
    }};
  }
};

//! An array containing pointers to all shape data types' angle function
constexpr auto angleFunctions = temple::tuples::unpackToFunction<
  allShapeDataTypes,
  angleFunctionFunctor
>();

} // namespace data

/* Core shape data, this has dynamic types and is hence initialized in the
 * .cpp file from the tuple containing all shape data types and the
 * shapeInformationFunctor
 */
const ShapeDataMapType& shapeData();


/*! Main data source of symmetries, derived from the various Symmetry data
 * classes. Use construct-on-first-use idiom to avoid static initialization
 * disasters accessing this variable
 *
 * Subtleties:
 * - By using an instance instead of a pointer-to-new, this is destructed
 *   properly before exit and does not "leak"
 * - This however introduces another possible issue if other static object
 *   destructors use this variable, because once again, the order of static
 *   deinitialization is random.
 */
const ShapeDataMapType& shapeData() {
  static const auto dataMap = temple::tuples::unpackToFunction<
    data::allShapeDataTypes,
    data::shapeInformationFunctor
  >();

  return dataMap;
}

AngleFunction angleFunction(const Shape shape) {
  auto shapeIndex = static_cast<unsigned>(shape);
  return [shapeIndex](const Vertex a, const Vertex b) {
    return data::angleFunctions.at(shapeIndex)(a, b);
  };
}

Coordinates coordinates(const Shape shape) {
  return shapeData().at(shape).coordinates;
}

std::string spaceFreeName(const Shape shape) {
  std::string toModify = shapeData().at(shape).stringName;

  std::replace(
    toModify.begin(),
    toModify.end(),
    ' ',
    '-'
  );

  return toModify;
}

unsigned nameIndex(const Shape shape) {
  return std::find(
    allShapes.begin(),
    allShapes.end(),
    shape
  ) - allShapes.begin();
}

const std::string& name(const Shape shape) {
  return shapeData().at(shape).stringName;
}

/*! @brief Fetch the shape name from its string
 *
 * @complexity{@math{\Theta(S)}}
 * @throws std::logic_error if no matching shape can be found
 */
Shape nameFromString(const std::string& shapeNameString) {
  for(const Shape shape : allShapes) {
    if(shapeData().at(shape).stringName == shapeNameString) {
      return shape;
    }
  }

  throw std::logic_error("No shape exists under that string name!");
}

unsigned size(const Shape shape) {
  return shapeData().at(shape).size;
}

const RotationsList& rotations(const Shape shape) {
  return shapeData().at(shape).rotations;
}

const Permutation& mirror(const Shape shape) {
  return shapeData().at(shape).mirror;
}

PointGroup pointGroup(const Shape shape) {
  return shapeData().at(shape).pointGroup;
}

const TetrahedronList& tetrahedra(const Shape shape) {
  return shapeData().at(shape).tetrahedra;
}

bool threeDimensional(const Shape shape) {
  return shapeData().at(shape).threeDimensional;
}

} // namespace shapes
} // namespace Scine
