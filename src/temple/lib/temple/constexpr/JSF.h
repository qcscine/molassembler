#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_JSF
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_JSF

#include <cstdint>
#include <random>

namespace temple {

namespace jsf {

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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

template<
  typename UnsignedType,
  unsigned p,
  unsigned q,
  unsigned r
> class JSF {
private:
  UnsignedType _a, _b, _c, _d;

  static constexpr unsigned bits = 8 * sizeof(UnsignedType);

  static constexpr UnsignedType _rotate(UnsignedType x, unsigned k) {
    return (x << k) | (x >> (bits - k));
  }

  constexpr void _advance() {
    UnsignedType e = _a - _rotate(_b, p);
    _a = _b ^ _rotate(_c, q);
    _b = _c + (r ? _rotate(_d, r) : _d);
    _c = _d + e;
    _d = e + _a;
  }

  constexpr void _advance(unsigned N) {
    for(unsigned i = 0; i < N; ++i) {
      _advance();
    }
  }

  constexpr void _seed(const std::array<UnsignedType, 4>& state) {
    _a = state[0];
    _b = state[1];
    _c = state[2];
    _d = state[3];
  }

  void _seed(std::seed_seq& seedSeq) {
    std::array<UnsignedType, 4> stateArray;

    seedSeq.generate(
      std::begin(stateArray),
      std::end(stateArray)
    );

    //! C++17 if constexpr
    if(std::is_same<UnsignedType, std::uint64_t>::value) {
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

    _seed(stateArray);
  }


public:
  using result_type = UnsignedType;

  static constexpr result_type min() {
    return 0;
  }

  static constexpr result_type max() {
    return ~result_type(0);
  }

  constexpr void seed(const std::array<UnsignedType, 4>& input) {
    _seed(input);
  }

  void seed(std::seed_seq& seedSeq) {
    _seed(seedSeq);
  }

  void seed(int seed) {
    std::seed_seq seedSeq {seed};
    _seed(seedSeq);
  }

  constexpr explicit JSF() = default;

  constexpr explicit JSF(const std::array<UnsignedType, 4>& input) {
    static_assert(
      std::is_unsigned<UnsignedType>::value,
      "The underlying type of the JSF generator must be unsigned!"
    );

    _seed(input);
  }

  explicit JSF(std::seed_seq& seedSeq) {
    _seed(seedSeq);
  }

  explicit JSF(int seed) {
    _seed(seed);
  }

  constexpr UnsignedType operator() () {
    _advance();

    return _d;
  }

  constexpr bool operator == (const JSF& other) const {
    return (
      _a == other._a
      && _b == other._b
      && _c == other._c
      && _d == other._d
    );
  }

  constexpr bool operator != (const JSF& other) const {
    return !(*this == other);
  }
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

} // namespace jsf

} // namespace temple

#endif // INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_JSF
