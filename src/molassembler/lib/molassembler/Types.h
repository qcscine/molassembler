// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_SHARED_TYPES_H
#define INCLUDE_MOLASSEMBLER_SHARED_TYPES_H

#include <cstddef>

/*!@file
 *
 * @brief Defines basic types widely shared across the project.
 */

namespace molassembler {

/*!
 * Bond type enumeration. Besides the classic organic single, double and triple
 * bonds, bond orders up to sextuple are explicitly included.
 */
enum class BondType : unsigned {
  Single,
  Double,
  Triple,
  Quadruple,
  Quintuple,
  Sextuple,
  Eta
};

enum class LengthUnit {
  Bohr,
  Angstrom
};

//! Unsigned integer atom index type. Used to refer to particular atoms.
using AtomIndex = std::size_t;

//! Type used to refer to particular bonds. Orders first < second.
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
