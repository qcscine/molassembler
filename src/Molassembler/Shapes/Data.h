/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Centralizes basic shape data in runtime types
 *
 * Central inclusion point of the library. Defines the main shape data and
 * all accessors. Shapes are defined in a constexpr fashion and homogenized
 * into a single container at compile-time, allowing both compile-time
 * computation and universal accessibility at run-time.
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_DATA_H
#define INCLUDE_MOLASSEMBLER_SHAPES_DATA_H

#include "boost/optional/optional_fwd.hpp"
#include "Eigen/Core"

#include "Molassembler/Shapes/Shapes.h"
#include "Molassembler/Shapes/PointGroups.h"
#include "Molassembler/Temple/StrongIndex.h"
#include "Molassembler/Temple/Preprocessor.h"

#include <vector>

namespace Scine {
namespace Molassembler {
namespace Shapes {

//! @brief Index of a shape vertex
struct vertex_index_tag;
using Vertex = Temple::StrongIndex<vertex_index_tag, unsigned>;

//! @brief Representation of a shape vertex permutation
using Permutation = std::vector<Vertex>;

/* Typedefs */
//! The type to store shape rotations
using RotationsList = std::vector<Permutation>;

/*! @brief Function returning angle between vertices
 *
 * All angle functions can be called with arbitrary (valid) parameters
 * without failing. Valid here means that a != b and less than the size of
 * the shape requested.
 *
 * They return angles in radians.
 */
using AngleFunction = std::function<
  double(const Vertex, const Vertex)
>;

/*! @brief  This is a four-vertex tetrahedron definition.
 *
 * All shapes use None to signal to replace this position with the central atom
 * as it is not part of the vertex indexing scheme.
 */
using Tetrahedron = std::array<boost::optional<Vertex>, 4>;

/*!
 * All shapes use boost::none to signal to replace this position with the
 * central atom as it is not part of the vertex indexing scheme.
 */
using TetrahedronList = std::vector<Tetrahedron>;

//! Representation of idealized shape coordinates (does not include centroid)
using Coordinates = Eigen::Matrix<double, 3, Eigen::Dynamic>;

namespace Detail {
template<size_t ... Inds>
constexpr auto makeAllShapes(std::index_sequence<Inds...> /* seq */) {
  return std::array<Shape, nShapes> {{ static_cast<Shape>(Inds)... }};
}
} // namespace Detail

//! A list of all the enum class values
constexpr std::array<Shape, nShapes> allShapes = Detail::makeAllShapes(
  std::make_index_sequence<nShapes>()
);

/* Interface */
/*! @brief Fetch the string name of a shape
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT const std::string& name(Shape shape);

/*! @brief Fetch the shape name from its string
 *
 * @complexity{@math{\Theta(S)}}
 * @throws std::logic_error if no matching shape can be found
 */
MASM_EXPORT Shape nameFromString(const std::string& shapeNameString);

/*! @brief Fetch a space-free string of a shape for file naming
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT std::string spaceFreeName(Shape shape);

/*! @brief Fetch the number of vertices of a shape
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT unsigned size(Shape shape);

/*! @brief Fetches a shape's list of rotations
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT const RotationsList& rotations(Shape shape);

/*! @brief Fetches the mirror index mapping for a particular shape
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT const Permutation& mirror(Shape shape);

/*! @brief Gets a shape's angle function
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT AngleFunction angleFunction(Shape shape);

/*! @brief Fetch a shape's idealized coordiantes
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT Coordinates coordinates(Shape shape);

/*! @brief Get a shape's point group
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT PointGroup pointGroup(Shape shape);

/*! @brief Returns the index of a shape within allShapes
 *
 * @complexity{@math{\Theta(S)}}
 */
MASM_EXPORT unsigned nameIndex(Shape shape);

/*! @brief Fetches the list of tetrahedra defined in a shape
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT const TetrahedronList& tetrahedra(Shape shape);

/*! @brief Returns whether a shape is three dimensional
 *
 * @complexity{@math{\Theta(1)}}
 */
MASM_EXPORT bool threeDimensional(Shape shape);

} // namespace Shapes
} // namespace Molassembler
} // namespace Scine

#endif
