/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Data struct for storing a numeric interval.
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_VALUE_BOUNDS_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_VALUE_BOUNDS_H

#include <stdexcept>

namespace Scine {

namespace molassembler {

namespace DistanceGeometry {

/* Exchanged types */
struct ValueBounds {
  double lower, upper;

  /**
   * @brief Main constructor
   *
   * @param passLower Lower end of the value boundaries
   * @param passUpper Upper end of the value boundaries
   *
   * @throws std::runtime_error If the upper value is smaller than the lower.
   *   Equalities are allowed.
   *
   * @return A ValueBounds instance with the specified bounds
   */
  constexpr ValueBounds(double passLower, double passUpper)
    : lower(passLower), upper(passUpper)
  {
    if(upper < lower) {
      throw std::runtime_error("Passed lower value is not smaller than the upper value!");
    }
  }

  /**
   * @brief Default constructor yields double lowest/max pair
   *
   * @return lowest/max bounds from numeric_limits
   */
  ValueBounds();

  /**
   * @brief Returns whether the bounds exactly equal another set of bounds
   *
   * @param other The other bounds to compare against
   *
   * @return If the bounds exactly equal one another
   */
  bool operator == (const ValueBounds& other) const;
  //! Negates @p operator ==
  bool operator != (const ValueBounds& other) const;
};

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine

#endif
