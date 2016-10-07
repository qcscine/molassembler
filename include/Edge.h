#ifndef INCLUDE_MOLECULE_EDGE_HPP
#define INCLUDE_MOLECULE_EDGE_HPP

#include <cassert>

#include "common_typedefs.h"

namespace MoleculeManip {

struct Edge {
  AtomIndexType i, j;
  BondType bondType;

  /* Constructors */
  Edge() = delete;
  Edge(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const BondType& bty
  ) : 
    i(std::min(a, b)), 
    j(std::max(a, b)), 
    bondType (bty) 
  {
    assert(a != b);
  };

  /* Information */
  bool smallerThan(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const {
    assert(a < b);
    return (
      i < a
      || (
        i == a
        && j < b 
      )
    );
  }

  bool greaterThan(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const {
    assert(a < b);
    return (
      i > a
      || (
        i == a
        && j > b 
      )
    );
  }

  /* Operators */
  bool operator < (const Edge& other) const {
    return (
      i < other.i
      || (
        i == other.i
        && j < other.j
      )
      || (
        i == other.i
        && j == other.j
        && bondType < other.bondType
      )
    );
  }
};

}

#endif
