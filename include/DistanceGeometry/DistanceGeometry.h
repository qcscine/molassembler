#ifndef DISTANCE_GEOMETRY_HPP
#define DISTANCE_GEOMETRY_HPP

#include <tuple>
#include "common_typedefs.h"
#include "DistanceGeometry/ValueBounds.h"

/* TODO
 */

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

/* Enum types */
/*! 
 * During distance geometry, when individual distance bounds are selected, 
 * other bounds must be adapted to the imposed restriction in metric space. The
 * set of algorithms to apply this adaptation are called metrization algorithms.
 */
enum class MetrizationOption {
  off,
//  partial, // -> unimplemented
  full
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
