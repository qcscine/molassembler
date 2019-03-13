/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Data struct to store chiral constraints for DG
 *
 * Contains some central data class declarations and type definitions for the
 * entire Distance Geometry scheme.
 */

#ifndef INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_H
#define INCLUDE_MOLASSEMBLER_DISTANCE_GEOMETRY_H

#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/Types.h"

#include <tuple>
#include <vector>
#include <cassert>
#include <array>

namespace Scine {

namespace molassembler {

//! Distance geometry-related classes and functions
namespace DistanceGeometry {

struct ChiralityConstraint {
  using AtomListType = std::vector<AtomIndex>;
  using SiteSequence = std::array<AtomListType, 4>;

  SiteSequence sites;
  double lower, upper;

  ChiralityConstraint(
    SiteSequence passSites,
    const double passLower,
    const double passUpper
  ) : sites(std::move(passSites)),
      lower(passLower),
      upper(passUpper)
  {
    // Must be <= because flat targets have lower = upper = 0
    assert(lower <= upper);
  }
};

struct DihedralConstraint {
  using AtomListType = std::vector<AtomIndex>;
  using SiteSequence = std::array<AtomListType, 4>;

  SiteSequence sites;
  double lower, upper;

  DihedralConstraint(
    SiteSequence passSites,
    const double passLower,
    const double passUpper
  ) : sites(std::move(passSites)),
      lower(passLower),
      upper(passUpper)
  {
    assert(lower <= upper);
  }
};

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine

#endif
