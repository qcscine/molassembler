/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "chemical_symmetries/Symmetries.h"

namespace Scine {

namespace Symmetry {

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
const SymmetryDataMapType& symmetryData() {
  static const auto dataMap = temple::TupleType::unpackToFunction<
    data::allShapeDataTypes,
    data::symmetryInformationFunctor
  >();

  return dataMap;
}

std::string spaceFreeName(const Shape shape) {
  std::string toModify = symmetryData().at(shape).stringName;

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

} // namespace Symmetry

} // namespace Scine
