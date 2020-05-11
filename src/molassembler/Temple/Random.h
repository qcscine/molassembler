/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Randomness helper functions
 *
 * Provides helpers for generating floating point, integer, boolean values.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_RANDOM_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_RANDOM_H

#include "molassembler/Temple/Traits.h"

#include <algorithm>
#include <array>
#include <random>

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Random {

//! Generate N floating point values in the range [lower, upper)
template<typename T, typename Engine>
std::enable_if_t<
  std::is_floating_point<T>::value,
  std::vector<T>
> getN(
  const T lower,
  const T upper,
  const unsigned N,
  Engine& engine
) {
  assert(lower <= upper);
  std::vector<T> returnNumbers;
  returnNumbers.reserve(N);

  std::uniform_real_distribution<T> uniformDistribution(lower, upper);

  for(unsigned i = 0; i < N; i++) {
    returnNumbers.emplace_back(uniformDistribution(engine));
  }

  return returnNumbers;
}

//! Generate N integral values in the range [lower, upper]
template<typename T, typename Engine>
std::enable_if_t<
  std::is_integral<T>::value,
  std::vector<T>
> getN(
  const T lower,
  const T upper,
  const unsigned N,
  Engine& engine
) {
  assert(lower <= upper);
  std::vector<T> returnNumbers;
  returnNumbers.reserve(N);

  std::uniform_int_distribution<T> uniformDistribution(lower, upper);

  for(unsigned i = 0; i < N; i++) {
    returnNumbers.emplace_back(uniformDistribution(engine));
  }

  return returnNumbers;
}

//! Generate a single floating point value in the range [lower, upper)
template<typename T, typename Engine>
std::enable_if_t<
  std::is_floating_point<T>::value,
  T
> getSingle(
  const T lower,
  const T upper,
  Engine& engine
) {
  assert(lower <= upper);
  std::uniform_real_distribution<T> uniformDistribution(lower, upper);
  return uniformDistribution(engine);
}

//! Generate a single integral value in the range [lower, upper]
template<typename T, typename Engine>
std::enable_if_t<
  std::is_integral<T>::value,
  T
> getSingle(
  const T lower,
  const T upper,
  Engine& engine
) {
  assert(lower <= upper);
  std::uniform_int_distribution<T> uniformDistribution(lower, upper);
  return uniformDistribution(engine);
}

//! Flip a coin
template<typename T, typename Engine>
std::enable_if_t<
  std::is_same<T, bool>::value,
  bool
> getSingle(Engine& engine) {
  std::uniform_int_distribution<unsigned> uniformDistribution(0, 1);
  return static_cast<bool>(
    uniformDistribution(engine)
  );
}

//! Pick a value using provided discrete distribution weights
template<class Container, typename Engine>
unsigned pickDiscrete(
  const Container& weights,
  const Engine& engine
) {
  std::discrete_distribution<unsigned> distribution {
    std::begin(weights),
    std::end(weights)
  };

  return distribution(engine);
}

/*! @brief Use underlying PRNG to shuffle a Container in-place
 *
 * @complexity{@math{\Theta(N)}}
 */
template<class Container, typename Engine>
void shuffle(
  Container& container,
  Engine& engine
) {
  std::shuffle(
    std::begin(container),
    std::end(container),
    engine
  );
}

} // namespace Random
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
