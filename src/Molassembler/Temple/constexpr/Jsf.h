/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Bob Jenkins' small fast noncryptographic PRNG
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_JSF
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_JSF

#include <cstdint>
#include <random>
#include <array>

namespace Scine {
namespace Molassembler {
namespace Temple {

/* Heavily modified from:
 *
 * A C++ implementation of a Bob Jenkins Small Fast (Noncryptographic) PRNGs
 * Based on code published by Bob Jenkins in 2007, adapted for C++
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2018 Melissa E. O'Neill
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/*!
 * @brief General class enabling the construction of a pattern of PRNGs by
 *   Bob Jenkins
 *
 * @tparam UnsignedType An unsigned integer type that contains the main state of
 *   the PRNG. Choose between 32 and 64-bit state here!
 * @tparam p A parameter that influences the state advance operations.
 * @tparam q A parameter that influences the state advance operations.
 * @tparam r A parameter that influences the state advance operations.
 *
 */
template<
  typename UnsignedType,
  unsigned p,
  unsigned q,
  unsigned r
> class JSF {
public:
  //! When used as a functor, this is the return type of operator ()
  using result_type = UnsignedType;

//!@name Public static properties
//!@{
  //! Minimum value of result_type
  static constexpr result_type min() {
    return 0;
  }

  //! Maximum value of result_type
  static constexpr result_type max() {
    return ~result_type(0);
  }
//!@}

//!@name Constructors
//!@{
  //! Default constructor
  constexpr explicit JSF() = default;

  //! Construct from four seed values
  constexpr explicit JSF(const std::array<UnsignedType, 4>& input) {
    static_assert(
      std::is_unsigned<UnsignedType>::value,
      "The underlying type of the JSF generator must be unsigned!"
    );

    seed_(input);
  }

  //! Construct from a seed sequence
  explicit JSF(std::seed_seq& seedSeq) {
    seed_(seedSeq);
  }

  //! Construct from a single integer seed value
  explicit JSF(int seedInt) {
    seed(seedInt);
  }
//!@}

//!@name Modifiers
//!@{
  //! Seed the underlying state with four values
  constexpr void seed(const std::array<UnsignedType, 4>& input) {
    seed_(input);
  }

  //! Seed the underlying state with a seed sequence
  void seed(std::seed_seq& seedSeq) {
    seed_(seedSeq);
  }

  //! Seed the underlying state with a single integer value
  void seed(int seed) {
    std::seed_seq seedSeq {{seed}};
    seed_(seedSeq);
  }
//!@}

//!@name Operators
//!@{
  /*! @brief Advance the state and return the current value
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr UnsignedType operator() () {
    advance_();

    return d_;
  }

  /*! @brief Compares the underlying state of two instances
   *
   * @complexity{@math{\Theta(1)}}
   */
  constexpr bool operator == (const JSF& other) const {
    return (
      a_ == other.a_
      && b_ == other.b_
      && c_ == other.c_
      && d_ == other.d_
    );
  }

  //! Compares the underlying state of two instances
  constexpr bool operator != (const JSF& other) const {
    return !(*this == other);
  }
//!@}

private:
//!@name State
//!@{
  UnsignedType a_, b_, c_, d_;
//!@}

//!@name Static properties
//!@{
  static constexpr unsigned bits = 8 * sizeof(UnsignedType);
//!@}

//!@name Private member functions
//!@{
  static constexpr UnsignedType rotate_(UnsignedType x, unsigned k) {
    return (x << k) | (x >> (bits - k));
  }

  constexpr void advance_() {
    UnsignedType e = a_ - rotate_(b_, p);
    a_ = b_ ^ rotate_(c_, q);
    b_ = c_ + ((r > 0) ? rotate_(d_, r) : d_);
    c_ = d_ + e;
    d_ = e + a_;
  }

  constexpr void advance_(unsigned N) {
    for(unsigned i = 0; i < N; ++i) {
      advance_();
    }
  }

  constexpr void seed_(const std::array<UnsignedType, 4>& state) {
    a_ = state[0];
    b_ = state[1];
    c_ = state[2];
    d_ = state[3];
  }

  void seed_(std::seed_seq& seedSeq) {
    std::array<UnsignedType, 4> stateArray;

    seedSeq.generate(
      std::begin(stateArray),
      std::end(stateArray)
    );

    /* You can safely ignore the warning that shift cout >= width of type for
     * instantiations with uint32_t by clang. In those cases, stateArray is
     * already adequately filled.
     */
    if constexpr(std::is_same<UnsignedType, std::uint64_t>::value) {
      /* seed_seq only generates 32 bit unsigneds, so just combine 8 32-bit
       * values into 4 64-bit values for the state array
       */
      std::array<UnsignedType, 4> topBits;

      seedSeq.generate(
        std::begin(topBits),
        std::end(topBits)
      );

      // Shift the topBits values by 32 and XOR them onto stateArray
      for(unsigned i = 0; i < 4; ++i) {
        stateArray[i] ^= (topBits[i] << 32);
      }
    }

    seed_(stateArray);
  }
//!@}
};

///// ---- Specific JSF Generators ---- ////
//
// Each size has variations corresponding to different parameter sets.
// Each variant will create a distinct (and hopefully statistically
// independent) sequence.
//
// The constants are all those suggested by Bob Jenkins.  The n variants
// perform only two rotations, the r variants perform three.
// 128 state bits, uint32_t output

using JSF32na = JSF<uint32_t, 27, 17, 0>;
using JSF32nb = JSF<uint32_t,  9, 16, 0>;
using JSF32nc = JSF<uint32_t,  9, 24, 0>;
using JSF32nd = JSF<uint32_t, 10, 16, 0>;
using JSF32ne = JSF<uint32_t, 10, 24, 0>;
using JSF32nf = JSF<uint32_t, 11, 16, 0>;
using JSF32ng = JSF<uint32_t, 11, 24, 0>;
using JSF32nh = JSF<uint32_t, 25,  8, 0>;
using JSF32ni = JSF<uint32_t, 25, 16, 0>;
using JSF32nj = JSF<uint32_t, 26,  8, 0>;
using JSF32nk = JSF<uint32_t, 26, 16, 0>;
using JSF32nl = JSF<uint32_t, 26, 17, 0>;
using JSF32nm = JSF<uint32_t, 27, 16, 0>;

using JSF32ra = JSF<uint32_t,  3, 14, 24>;
using JSF32rb = JSF<uint32_t,  3, 25, 15>;
using JSF32rc = JSF<uint32_t,  4, 15, 24>;
using JSF32rd = JSF<uint32_t,  6, 16, 28>;
using JSF32re = JSF<uint32_t,  7, 16, 27>;
using JSF32rf = JSF<uint32_t,  8, 14,  3>;
using JSF32rg = JSF<uint32_t, 11, 16, 23>;
using JSF32rh = JSF<uint32_t, 12, 16, 22>;
using JSF32ri = JSF<uint32_t, 12, 17, 23>;
using JSF32rj = JSF<uint32_t, 13, 16, 22>;
using JSF32rk = JSF<uint32_t, 15, 25,  3>;
using JSF32rl = JSF<uint32_t, 16,  9,  3>;
using JSF32rm = JSF<uint32_t, 17,  9,  3>;
using JSF32rn = JSF<uint32_t, 17, 27,  7>;
using JSF32ro = JSF<uint32_t, 19,  7,  3>;
using JSF32rp = JSF<uint32_t, 23, 15, 11>;
using JSF32rq = JSF<uint32_t, 23, 16, 11>;
using JSF32rr = JSF<uint32_t, 23, 17, 11>;
using JSF32rs = JSF<uint32_t, 24,  3, 16>;
using JSF32rt = JSF<uint32_t, 24,  4, 16>;
using JSF32ru = JSF<uint32_t, 25, 14,  3>;
using JSF32rv = JSF<uint32_t, 27, 16,  6>;
using JSF32rw = JSF<uint32_t, 27, 16,  7>;

using JSF32n = JSF32na;
using JSF32r = JSF32rq;
using JSF32  = JSF32n;

// - 256 state bits, uint64_t output

using JSF64na = JSF<uint64_t, 39, 11,  0>;
using JSF64ra = JSF<uint64_t,  7, 13, 37>;

using JSF64n = JSF64na;
using JSF64r = JSF64ra;
using JSF64  = JSF64r;

//! Provides a seeded engine instance
template<class Engine = JSF32>
struct Generator {
  Engine engine;

  explicit Generator() {
#ifndef NDEBUG
    engine.seed(272181374);
#else
    std::random_device randomDevice;

    engine.seed(
      std::array<typename Engine::result_type, 4> {
        randomDevice(),
        randomDevice(),
        randomDevice(),
        randomDevice()
      }
    );
#endif
  }
};

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
