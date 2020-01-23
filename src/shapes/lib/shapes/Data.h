/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Centralizes basic shape data in runtime types
 *
 * Central inclusion point of the library. Defines the main shape data and
 * all accessors. Shapes are defined in a constexpr fashion and homogenized
 * into a single container at compile-time, allowing both compile-time
 * computation and universal accessibility at run-time.
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_DATA_H
#define INCLUDE_MOLASSEMBLER_SHAPES_DATA_H

#include "boost/optional.hpp"
#include "Eigen/Core"
#include "temple/constexpr/TupleType.h"

#include "shapes/constexpr/Data.h"

#include <vector>
#include <unordered_map>

namespace Scine {

namespace shapes {

/* Typedefs */
//! The type to store shape rotations
using RotationsList = std::vector<
  std::vector<unsigned>
>;

/*!
 * All angle functions can be called with arbitrary (valid) parameters
 * without failing. Valid here means that a != b and less than the size of
 * the shape requested.
 *
 * They return angles in radians.
 */
using AngleFunctionType = std::function<
  double(const unsigned, const unsigned)
>;

/*!
 * All symmetries have a guess implementation of what could work as the defined
 * tetrahedra. Have to use boost::none to signal to replace this position with
 * the central atom as it is not part of the indexing scheme used here.
 *
 * In case all higher symmetries than trigonal pyramidal are representable
 * without boost::none and that proves to work, then perhaps make an exception
 * for it and treat all others without the optional. If that cannot be done,
 * consider refactoring (changing the numbering scheme in some fashion that
 * boost::none does not have to be used.
 */
using TetrahedronList = std::vector<
  std::array<
    boost::optional<unsigned>,
    4
  >
>;

using CoordinateList = Eigen::Matrix<double, 3, Eigen::Dynamic>;
using MirrorMap = std::vector<unsigned>;

//! Dynamic shape information data struct
struct ShapeInformation {
  const std::string stringName;
  const unsigned size;
  const RotationsList rotations;
  const TetrahedronList tetrahedra;
  const CoordinateList coordinates;
  const MirrorMap mirror;
  const PointGroup pointGroup;
  bool threeDimensional;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// Helper function to create all names vector
template<size_t ... Inds>
constexpr std::array<Shape, nShapes> makeAllShapes(
  std::index_sequence<Inds...> /* indexSequence */
) {
  return {{
    static_cast<Shape>(Inds)...
  }};
}

//! A list of all the enum class values
constexpr std::array<Shape, nShapes> allShapes = makeAllShapes(
  std::make_index_sequence<nShapes>()
);

//! Map type used to store shape information structs
using ShapeDataMapType = std::unordered_map<
  Shape,
  ShapeInformation,
  std::hash<Shape>,
  std::equal_to<Shape>,
  Eigen::aligned_allocator<std::pair<const Shape, ShapeInformation>>
>;

namespace data {

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
        tetrahedra.at(i).at(j) = constexprValue;
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
CoordinateList makeCoordinates(
  const std::array<temple::Vector, shapeSize>& constexprCoordinates
) {
  CoordinateList coordinates(3, shapeSize);

  for(unsigned i = 0; i < shapeSize; ++i) {
    coordinates.col(i) = toEigen(constexprCoordinates.at(i));
  }

  return coordinates;
}

template<size_t shapeSize>
MirrorMap makeMirror(
  const std::array<unsigned, shapeSize>& constexprMirror
) {
  std::vector<unsigned> mirror (shapeSize);
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

/* Interface */
/*! @brief Fetch the string name of a shape
 *
 * @complexity{@math{\Theta(1)}}
 */
inline const std::string& name(const Shape shape) {
  return shapeData().at(shape).stringName;
}

/*! @brief Fetch the shape name from its string
 *
 * @complexity{@math{\Theta(S)}}
 * @throws std::logic_error if no matching shape can be found
 */
inline Shape nameFromString(const std::string& shapeNameString) {
  for(const Shape shape : allShapes) {
    if(shapeData().at(shape).stringName == shapeNameString) {
      return shape;
    }
  }

  throw std::logic_error("No shape exists under that string name!");
}

/*! @brief Fetch a space-free string of a shape for file naming
 *
 * @complexity{@math{\Theta(1)}}
 */
std::string spaceFreeName(Shape shape);

/*! @brief Fetch the number of vertices of a shape
 *
 * @complexity{@math{\Theta(1)}}
 */
inline unsigned size(const Shape shape) {
  return shapeData().at(shape).size;
}

/*! @brief Fetches a shape's list of rotations
 *
 * @complexity{@math{\Theta(1)}}
 */
inline const RotationsList& rotations(const Shape shape) {
  return shapeData().at(shape).rotations;
}

/*! @brief Fetches the mirror index mapping for a particular shape
 *
 * @complexity{@math{\Theta(1)}}
 */
inline const MirrorMap& mirror(const Shape shape) {
  return shapeData().at(shape).mirror;
}

/*! @brief Gets a shape's angle function
 *
 * @complexity{@math{\Theta(1)}}
 */
inline data::AngleFunctionPtr angleFunction(const Shape shape) {
  auto shapeIndex = static_cast<unsigned>(shape);
  return data::angleFunctions.at(shapeIndex);
}

/*! @brief Get a shape's point group
 *
 * @complexity{@math{\Theta(1)}}
 */
inline PointGroup pointGroup(const Shape shape) {
  return shapeData().at(shape).pointGroup;
}

/*! @brief Returns the index of a shape within allShapes
 *
 * @complexity{@math{\Theta(S)}}
 */
PURITY_STRONG unsigned nameIndex(Shape shape);

/*! @brief Fetches the list of tetrahedra defined in a shape
 *
 * @complexity{@math{\Theta(1)}}
 */
inline const TetrahedronList& tetrahedra(const Shape shape) {
  return shapeData().at(shape).tetrahedra;
}

/*! @brief Returns whether a shape is three dimensional
 *
 * @complexity{@math{\Theta(1)}}
 */
inline bool threeDimensional(const Shape shape) {
  return shapeData().at(shape).threeDimensional;
}

} // namespace shapes

} // namespace Scine

#endif
