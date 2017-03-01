#ifndef INCLUDE_COMMON_TYPEDEFS_H
#define INCLUDE_COMMON_TYPEDEFS_H

#include <tuple>

namespace MoleculeManip {

/* Common typedefs */
enum class BondType : unsigned {
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Sextuple,
    Aromatic,
    Eta
};

using unsigned_type = unsigned;

using EdgeIndexType = unsigned;
using AtomIndexType = unsigned; 

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
