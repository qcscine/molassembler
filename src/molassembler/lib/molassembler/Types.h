#ifndef INCLUDE_MOLASSEMBLER_SHARED_TYPES_H
#define INCLUDE_MOLASSEMBLER_SHARED_TYPES_H

#include <cstddef>

namespace molassembler {

/*!
 * Bond type enumeration. Besides the classic organic single, double and triple
 * bonds, bond orders up to sextuple are explicitly included.
 *
 * Although currently unused, Aromatic and Eta bonds are included in
 * anticipation of their necessity.
 */
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

enum class LengthUnit {
  Bohr,
  Angstrom
};

using AtomIndex = std::size_t;

struct BondIndex {
  AtomIndex first, second;

  BondIndex();
  BondIndex(AtomIndex a, AtomIndex b) noexcept;

  bool operator < (const BondIndex& other) const;
  bool operator == (const BondIndex& other) const;
};

//! Descriptive name for dlib indices
using dlibIndexType = long;

/*! For bitmasks grouping components of immediate atom environments
 *
 * Differing strictnesses of comparisons may be desirable for various
 * purposes, hence a modular comparison function is provided.
 */
enum class AtomEnvironmentComponents : unsigned {
  ElementTypes,
  BondOrders,
  Symmetries,
  Stereopermutations // Symmetries must be set in conjunction with this
};

} // namespace molassembler

#endif
