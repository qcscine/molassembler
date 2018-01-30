#ifndef DISTANCE_GEOMETRY_HPP
#define DISTANCE_GEOMETRY_HPP

#include <tuple>
#include "common_typedefs.h"
#include "DistanceGeometry/ValueBounds.h"

/*! @file
 *
 * Contains some central data class declarations and type definitions for the
 * entire Distance Geometry scheme.
 */

namespace MoleculeManip {

namespace DistanceGeometry {

struct ChiralityConstraint {
  std::array<AtomIndexType, 4> indices;
  double lower, upper;

  ChiralityConstraint(
    const std::array<AtomIndexType, 4>& indices,
    const double& lower,
    const double& upper
  );
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
