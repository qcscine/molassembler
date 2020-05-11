/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Fixed-size bitset class
 *
 * Contains a fixed-size bitset implementation.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_BITSET_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_BITSET_H

#include "molassembler/Temple/constexpr/Math.h"
#include "molassembler/Temple/constexpr/Array.h"

namespace Scine {
namespace Temple {

/**
 * @brief Constexpr bitset class
 *
 * @tparam N Number of bits to store.
 */
template<std::size_t N>
class Bitset {
public:
//!@name Constructor
//!@{
  /*! @brief Zeroing constructor
   *
   * @complexity{@math{\Theta(N)}}
   */
  explicit constexpr Bitset() { zero(); }
//!@}

//!@name Modification
//!@{
  /*! @brief Zero out all bits
   *
   * @complexity{@math{\Theta(N)}}
   */
  constexpr void zero() {
    for(auto& block : storage) {
      block = 0;
    }
  }

  /*! @brief Sets a specific bit
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr void set(const std::size_t i) {
    std::size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    std::size_t bitIndex = i - bitsPerBlock * blockIndex;

    storage.at(blockIndex) |= (1ull << bitIndex);
  }

  /*! @brief Unsets a specific bit
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr void unset(const std::size_t i) {
    std::size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    std::size_t bitIndex = i - bitsPerBlock * blockIndex;

    storage.at(blockIndex) ^= (1ull << bitIndex);
  }

  /*! @brief Sets a specific bit to a specified value
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr void set(const std::size_t i, const bool value) {
    if(value) {
      set(i);
    } else {
      unset(i);
    }
  }
//!@}

//!@name Information
//!@{
  /*! @brief Test the value at a particular bit
   *
   * @complexity{@math{\Theta(1)}}
   */
  PURITY_STRONG constexpr bool test(const std::size_t i) const {
    std::size_t blockIndex = Math::floor(static_cast<double>(i) / bitsPerBlock);
    std::size_t bitIndex = i - bitsPerBlock * blockIndex;

    return (
      storage.at(blockIndex) & (1ull << bitIndex)
    ) != 0;
  }
//!@}

private:
//!@name Types
//!@{
  using Block = long long unsigned;
//!@}

//!@name Static const values
//!@{
  //! Number of bits per block
  static constexpr std::size_t bitsPerBlock = Math::floor(
    Math::log(
      static_cast<double>(std::numeric_limits<Block>::max()),
      2.0
    )
  );

  //! Number of Block types stored
  static constexpr std::size_t B = Math::ceil(
    static_cast<double>(N) / static_cast<double>(bitsPerBlock)
  );
//!@}

//!@name State
//!@{
  Array<Block, B> storage;
//!@}

};

} // namespace Temple
} // namespace Scine

#endif
