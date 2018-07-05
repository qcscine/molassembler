#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_BITSET_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_BITSET_H

#include "Math.h"
#include "Array.h"

/*!@file
 *
 * Contains a fixed-size bitset implementation.
 */

namespace temple {

template<size_t N>
struct Bitset {
  using Block = long long unsigned;

  static constexpr size_t bitsPerBlock = Math::floor(
    Math::log(
      static_cast<double>(std::numeric_limits<Block>::max()),
      2.0
    )
  );

  static constexpr size_t B = Math::ceil(
    static_cast<double>(N) / static_cast<double>(bitsPerBlock)
  );

  Array<Block, B> storage;

  explicit constexpr Bitset() { zero(); }

  //! Zero out all bits
  constexpr void zero() {
    for(auto& block : storage) {
      block = 0;
    }
  }

  //! Sets a specific bit
  constexpr void set(size_t i) {
    size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    size_t bitIndex = i - bitsPerBlock * blockIndex;

    storage.at(blockIndex) |= (1ull << bitIndex);
  }

  //! Unsets a specific bit
  constexpr void unset(size_t i) {
    size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    size_t bitIndex = i - bitsPerBlock * blockIndex;

    storage.at(blockIndex) ^= (1ull << bitIndex);
  }

  //! Sets a specific bit to a specified value
  constexpr void set(size_t i, bool value) {
    if(value) {
      set(i);
    } else {
      unset(i);
    }
  }

  constexpr bool test(size_t i) const PURITY_STRONG {
    size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    size_t bitIndex = i - bitsPerBlock * blockIndex;

    return (
      storage.at(blockIndex) & (1ull << bitIndex)
    ) != 0;
  }
};

} // namespace temple

#endif
