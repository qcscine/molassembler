#ifndef INCLUDE_COMMON_TYPEDEFS_H
#define INCLUDE_COMMON_TYPEDEFS_H

#include <tuple>

namespace MoleculeManip {

/* Common typedefs */
enum class BondType : uint8_t {
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Sextuple,
    Aromatic
};

using unsigned_type = unsigned;
using EdgeIndexType = uint32_t;
using AtomIndexType = uint16_t; // 65k max is sufficient

/* derived types */
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

} // eo namespace

#endif
