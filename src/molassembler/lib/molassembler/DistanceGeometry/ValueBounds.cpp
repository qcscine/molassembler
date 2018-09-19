// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/DistanceGeometry/ValueBounds.h"

#include <limits>
#include <cassert>
#include <utility>

namespace molassembler {

namespace DistanceGeometry {

ValueBounds::ValueBounds() : ValueBounds(
  std::numeric_limits<double>::lowest(),
  std::numeric_limits<double>::max()
) {}

ValueBounds::ValueBounds(
  const double passLower,
  const double passUpper
) : lower(passLower),
    upper(passUpper)
{
  assert(lower <= upper);
}

} // namespace molassembler

} // namespace DistanceGeometry
