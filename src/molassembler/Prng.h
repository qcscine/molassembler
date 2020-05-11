/*!@file
 * @brief Engine wrapper around temple's JSF PRNG for centralized re-seeding
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_PRNG_H
#define INCLUDE_MOLASSEMBLER_PRNG_H

#include "molassembler/Export.h"
#include <memory>
#include <vector>

namespace Scine {
namespace Molassembler {

/**
 * @brief Randomness source for the library
 */
namespace Random {

//! @brief Drives a PRNG
class MASM_EXPORT Engine {
public:
  //! The type this engine generates
  using result_type = uint32_t;

  explicit Engine();
  explicit Engine(int seedArg);
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

  /*! @brief Advances the state and returns a value
   *
   * @complexity{@math{\Theta(1)}}
   */
  result_type operator() () const;

  //! Compare this engine's state with that of another engine
  bool operator == (const Engine& other) const;

private:
  struct Impl;
  std::unique_ptr<Impl> pImpl_;
};

} // namespace Random
} // namespace Molassembler
} // namespace Scine

#endif
