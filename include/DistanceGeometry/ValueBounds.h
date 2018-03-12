#ifndef INCLUDE_DISTANCE_GEOMETRY_VALUE_BOUNDS_H
#define INCLUDE_DISTANCE_GEOMETRY_VALUE_BOUNDS_H

/*! @file
 *
 * Defines a data struct for storing a numeric interval.
 */

namespace molassembler {

namespace DistanceGeometry {

/* Exchanged types */
struct ValueBounds {
  double lower, upper;

  ValueBounds();

  ValueBounds(const double& lower, const double& upper);

  ValueBounds(const ValueBounds& other);

  ValueBounds& operator = (const ValueBounds& other);
  ValueBounds& operator = (ValueBounds&& other);
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
