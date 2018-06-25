#include "DistanceGeometry/ValueBounds.h"

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
  const double lower,
  const double upper
) : lower(lower),
    upper(upper)
{
  assert(lower <= upper);
}

ValueBounds::ValueBounds(const ValueBounds& other) : ValueBounds(other.lower, other.upper) {}
ValueBounds& ValueBounds::operator = (const ValueBounds& other) = default;
ValueBounds& ValueBounds::operator = (ValueBounds&& other) noexcept = default;

} // namespace molassembler

} // namespace DistanceGeometry
