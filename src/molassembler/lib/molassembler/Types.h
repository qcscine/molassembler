/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Defines basic types widely shared across the project.
 */

#ifndef INCLUDE_MOLASSEMBLER_SHARED_TYPES_H
#define INCLUDE_MOLASSEMBLER_SHARED_TYPES_H

#include <cstddef>
#include <iterator>

namespace Scine {

namespace molassembler {

/*!
 * @brief Discrete bond type numeration
 *
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
  /*! @brief Internal bond order to mark haptic binding sites
   *
   * Bond order used internally only that relabels bonds in haptic binding
   * sites. Name from the standard eta notation for haptic ligands.
   */
  Eta
};

constexpr unsigned nBondTypes = 7;

enum class LengthUnit {
  Bohr,
  Angstrom
};

//! Unsigned integer atom index type. Used to refer to particular atoms.
using AtomIndex = std::size_t;

//! Type used to refer to particular bonds. Orders first < second.
struct BondIndex {
  /*! Iterator type is just a pointer to member
   *
   * Since the standard guarantees that struct members are laid out
   * successively in memory, we can implement an iterator using just pointers.
   * It doesn't have the necessary typedefs for STL algorithms, but could
   * be useful just the same.
   */
  using const_iterator = const AtomIndex*;

  AtomIndex first, second;

  BondIndex();
  BondIndex(AtomIndex a, AtomIndex b) noexcept;

  //! Whether first or second is @p a
  bool contains(AtomIndex a) const;

  //! Lexicographic comparison
  bool operator < (const BondIndex& other) const;
  //! Lexicographic comparison
  bool operator == (const BondIndex& other) const;

  const_iterator begin() const {
    return &first;
  }

  const_iterator end() const {
    return std::next(&second);
  }
};

/*! @brief Hash for BondIndex so it can be used as a key type in unordered containers
 *
 * @complexity{@math{\Theta(1)}}
 */
std::size_t hash_value(const BondIndex& bond);

/*!
 * @brief For bitmasks grouping components of immediate atom environments
 *
 * Differing strictnesses of comparisons may be desirable for various
 * purposes, hence a modular comparison function is provided.
 *
 * @warning Setting Stereopermutations without setting Shapes does nothing.
 */
enum class AtomEnvironmentComponents : unsigned {
  None = 0,
  ElementTypes = (1 << 0),
  BondOrders = (1 << 1),
  Shapes = (1 << 2),
  Stereopermutations = (1 << 3),
  All = ElementTypes | BondOrders | Shapes | Stereopermutations
};

} // namespace molassembler

} // namespace Scine

/* Operators for bitmask-like manipulation of AtomEnvironmentComponents must be
 * at global scope, otherwise they can interfere with name lookup.
 */
constexpr inline bool operator & (
  Scine::molassembler::AtomEnvironmentComponents a,
  Scine::molassembler::AtomEnvironmentComponents b
) {
  using T = Scine::molassembler::AtomEnvironmentComponents;
  return (
    static_cast<std::underlying_type_t<T>>(a)
    & static_cast<std::underlying_type_t<T>>(b)
  ) != 0;
}

constexpr inline Scine::molassembler::AtomEnvironmentComponents operator | (
  Scine::molassembler::AtomEnvironmentComponents a,
  Scine::molassembler::AtomEnvironmentComponents b
) {
  using T = Scine::molassembler::AtomEnvironmentComponents;

  return static_cast<T>(
    static_cast<std::underlying_type_t<T>>(a)
    | static_cast<std::underlying_type_t<T>>(b)
  );
}

#endif
