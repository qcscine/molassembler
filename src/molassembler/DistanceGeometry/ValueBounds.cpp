/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/DistanceGeometry/ValueBounds.h"

#include <limits>

namespace Scine {

namespace Molassembler {

namespace DistanceGeometry {

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

} // namespace Molassembler

} // namespace DistanceGeometry

} // namespace Scine
