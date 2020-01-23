/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "shapes/Data.h"

namespace Scine {
namespace shapes {
namespace data {

Eigen::Vector3d toEigen(const temple::Vector& cVector) {
  return {
    cVector.data[0],
    cVector.data[1],
    cVector.data[2]
  };
}

} // namespace data

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

} // namespace shapes
} // namespace Scine
