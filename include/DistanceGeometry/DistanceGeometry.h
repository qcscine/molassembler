#ifndef DISTANCE_GEOMETRY_HPP
#define DISTANCE_GEOMETRY_HPP

#include <tuple>
#include "AdjacencyList.h"

/* TODO
 */

namespace MoleculeManip {

namespace DistanceGeometry {

/* Typedefs */
struct ChiralityConstraint {
  std::array<AtomIndexType, 4> indices;
  double lower, upper;

  ChiralityConstraint(
    const std::array<AtomIndexType, 4>& indices,
    const double& lower,
    const double& upper
  ) : indices(indices),
      lower(lower),
      upper(upper)
  {
    assert(lower < upper);
  }
};

/* Enum types */
enum class MetrizationOption {
  off,
//  partial, // -> unimplemented
  full
};

} // namespace DistanceGeometry

} // namespace MoleculeManip

#endif
