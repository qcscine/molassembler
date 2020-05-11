/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Data struct for storing a numeric interval.
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_VALUE_BOUNDS_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_VALUE_BOUNDS_H

#include <stdexcept>

namespace Scine {

namespace Molassembler {

namespace DistanceGeometry {

/* Exchanged types */
/**
 * @brief Data class for bounded values
 */
struct ValueBounds {
  double lower, upper;

  /**
   * @brief Main constructor
   *
   * @complexity{@math{\Theta(1)}}
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

  /** @brief Default constructor yields double lowest/max pair
   *
   * @complexity{@math{\Theta(1)}}
   * @return lowest/max bounds from numeric_limits
   */
  ValueBounds();

  /** @brief Returns whether the bounds exactly equal another set of bounds
   *
   * @complexity{@math{\Theta(1)}}
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

} // namespace Molassembler

} // namespace Scine

#endif
