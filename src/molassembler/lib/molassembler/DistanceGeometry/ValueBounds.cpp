/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/DistanceGeometry/ValueBounds.h"

#include <limits>

namespace Scine {

namespace molassembler {

namespace distance_geometry {

ValueBounds::ValueBounds() : ValueBounds(
  std::numeric_limits<double>::lowest(),
  std::numeric_limits<double>::max()
) {}

bool ValueBounds::operator == (const ValueBounds& other) const {
  return (
    lower == other.lower
    && upper == other.upper
  );
}

bool ValueBounds::operator != (const ValueBounds& other) const {
  return !(*this == other);
}

} // namespace molassembler

} // namespace distance_geometry

} // namespace Scine
