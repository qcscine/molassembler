#include "Symmetries.h"

namespace Symmetry {

namespace data {

Eigen::Vector3d toEigen(const constable::Vector& cVector) {
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
const std::map<Name, SymmetryInformation>& symmetryData() {
  static const auto dataMap = constable::TupleType::unpackToFunction<
    data::allSymmetryDataTypes,
    data::symmetryInformationFunctor
  >();

  return dataMap;
}

} // namespace Symmetry
