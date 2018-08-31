#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_VALUE_BOUNDS_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_VALUE_BOUNDS_H

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
  ValueBounds(double passLower, double passUpper);
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
