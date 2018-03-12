#include "DistanceGeometry/DistanceGeometry.h"

namespace molassembler {

namespace DistanceGeometry {

ValueBounds::ValueBounds() : ValueBounds(
  std::numeric_limits<double>::lowest(),
  std::numeric_limits<double>::max()
) {}

ValueBounds::ValueBounds(
  const double& lower,
  const double& upper
) : lower(lower),
    upper(upper)
{
  assert(lower <= upper);
}

ValueBounds::ValueBounds(const ValueBounds& other) : ValueBounds(other.lower, other.upper) {}

ValueBounds& ValueBounds::operator = (const ValueBounds& other) {
  lower = other.lower;
  upper = other.upper;
  return *this;
}

ValueBounds& ValueBounds::operator = (ValueBounds&& other) {
  std::swap(lower, other.lower);
  std::swap(upper, other.upper);
  return *this;
}


ChiralityConstraint::ChiralityConstraint(
  const std::array<AtomIndexType, 4>& indices,
  const double& lower,
  const double& upper
) : indices(indices),
    lower(lower),
    upper(upper)
{
  // Must be <= because flat targets have lower = upper = 0
  assert(lower <= upper);
}

} // namespace molassembler

} // namespace DistanceGeometry

