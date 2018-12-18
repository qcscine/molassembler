/*!@file
 * @brief Engine wrapper around temple's JSF PRNG for centralized re-seeding
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#ifndef INCLUDE_MOLASSEMBLER_PRNG_H
#define INCLUDE_MOLASSEMBLER_PRNG_H

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

#include <memory>
#include <vector>
#include <array>
#include <algorithm>
#include <random>

namespace Scine {

namespace molassembler {

/**
 * @brief Randomness source for the library
 */
namespace random {

//! Instances of this class can be used to drive PRNGs
class Engine {
public:
  //! The type this engine generates
  using result_type = uint32_t;

  explicit Engine();
  Engine(Engine&& other) noexcept;
  Engine& operator = (Engine&& other) noexcept;
  Engine(const Engine& other);
  Engine& operator = (const Engine& other);
  ~Engine();

  //! Minimum value of result_type
  static constexpr result_type min() {
    return 0;
  }

  //! Maximum value of result_type
  static constexpr result_type max() {
    return ~result_type(0);
  }

  //! Seed the underlying state with an integer value
  void seed(int x);

  //! Seed the underlying state with multiple integer values
  void seed(const std::vector<int>& signedSeeds);

  //! Advances the state and returns a value
  result_type operator() () const;

private:
  struct Impl;

#ifdef MOLASSEMBLER_ENABLE_PROPAGATE_CONST
  std::experimental::propagate_const<
    std::unique_ptr<Impl>
  > _pImpl;
#else
  std::unique_ptr<Impl> _pImpl;
#endif
};

} // namespace random

} // namespace molassembler

} // namespace Scine

#endif
