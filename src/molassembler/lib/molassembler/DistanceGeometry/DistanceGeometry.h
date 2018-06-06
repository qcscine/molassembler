#ifndef DISTANCE_GEOMETRY_HPP
#define DISTANCE_GEOMETRY_HPP

#include <tuple>
#include "detail/SharedTypes.h"
#include "DistanceGeometry/ValueBounds.h"

/*! @file
 *
 * Contains some central data class declarations and type definitions for the
 * entire Distance Geometry scheme.
 */

namespace molassembler {

//! Distance geometry-related classes and functions
namespace DistanceGeometry {

struct ChiralityConstraint {
  using AtomListType = std::vector<AtomIndexType>;
  using LigandSequence = std::array<AtomListType, 4>;

  LigandSequence sites;
  double lower, upper;

  ChiralityConstraint(
    const LigandSequence& sites,
    const double lower,
    const double upper
  ) : sites(sites),
      lower(lower),
      upper(upper)
  {
    // Must be <= because flat targets have lower = upper = 0
    assert(lower <= upper);
  }
};

enum class Partiality {
  FourAtom,
  TenPercent,
  All
};

} // namespace DistanceGeometry

} // namespace molassembler

#endif
