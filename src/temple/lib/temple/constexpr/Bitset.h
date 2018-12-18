/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Fixed-size bitset class
 *
 * Contains a fixed-size bitset implementation.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_BITSET_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_BITSET_H

#include "temple/constexpr/Math.h"
#include "temple/constexpr/Array.h"

namespace temple {

template<std::size_t N>
struct Bitset {
//!@name Types
//!@{
  using Block = long long unsigned;
//!@}

//!@name Static const values
//!@{
  static constexpr std::size_t bitsPerBlock = Math::floor(
    Math::log(
      static_cast<double>(std::numeric_limits<Block>::max()),
      2.0
    )
  );

  static constexpr std::size_t B = Math::ceil(
    static_cast<double>(N) / static_cast<double>(bitsPerBlock)
  );
//!@}

//!@name State
//!@{
  Array<Block, B> storage;
//!@}

//!@name Constructor
//!@{
  explicit constexpr Bitset() { zero(); }
//!@}

//!@name Modification
//!@{
  //! Zero out all bits
  constexpr void zero() {
    for(auto& block : storage) {
      block = 0;
    }
  }

  //! Sets a specific bit
  constexpr void set(std::size_t i) {
    std::size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    std::size_t bitIndex = i - bitsPerBlock * blockIndex;

    storage.at(blockIndex) |= (1ull << bitIndex);
  }

  //! Unsets a specific bit
  constexpr void unset(std::size_t i) {
    std::size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    std::size_t bitIndex = i - bitsPerBlock * blockIndex;

    storage.at(blockIndex) ^= (1ull << bitIndex);
  }

  //! Sets a specific bit to a specified value
  constexpr void set(std::size_t i, bool value) {
    if(value) {
      set(i);
    } else {
      unset(i);
    }
  }
//!@}

//!@name Information
//!@{
  PURITY_STRONG constexpr bool test(std::size_t i) const {
    std::size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    std::size_t bitIndex = i - bitsPerBlock * blockIndex;

    return (
      storage.at(blockIndex) & (1ull << bitIndex)
    ) != 0;
  }
//!@}
};

} // namespace temple

#endif
