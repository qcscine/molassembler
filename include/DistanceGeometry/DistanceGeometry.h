#ifndef DISTANCE_GEOMETRY_HPP
#define DISTANCE_GEOMETRY_HPP

#include <tuple>
#include "AdjacencyList.h"

/* TODO
 * - alter DistanceConstraint and ChiralityConstraint to structs
 */

namespace MoleculeManip {

namespace DistanceGeometry {

/* Typedefs */
using DistanceConstraint = std::tuple<
  AtomIndexType, // i
  AtomIndexType, // j
  double, // lower
  double // upper
>;

using ChiralityConstraint = std::tuple<
  AtomIndexType, // i
  AtomIndexType, // j
  AtomIndexType, // k
  AtomIndexType, // l
  double // target
>;

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
