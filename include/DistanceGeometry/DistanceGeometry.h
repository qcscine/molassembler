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
  double target;

  ChiralityConstraint(
    std::array<AtomIndexType, 4>&& indices,
    double&& target
  ) : indices(indices),
      target(target)
  {}

  ChiralityConstraint(
    const std::array<AtomIndexType, 4>& indices,
    const double& target
  ) : indices(indices),
      target(target)
  {}
};

/* Enum types */
enum class MetrizationOption {
  off,
//  partial, // -> unimplemented
  full
};

enum class EmbeddingOption : unsigned {
  threeDimensional = 3,
  fourDimensional = 4
};

} // eo namespace DistanceGeometry

} // eo namespace MoleculeManip

#endif
