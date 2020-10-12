/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Prng.h"

#include "Molassembler/Temple/constexpr/Jsf.h"
#include <random>

namespace Scine {
namespace Molassembler {
namespace Random {

struct Engine::Impl {
  Temple::JSF32 engine;

  void seed(int x) {
    engine.seed(x);
  }

  void seed(const std::vector<int>& seeds) {
    std::seed_seq seedSeq(
      std::begin(seeds),
      std::end(seeds)
    );

    engine.seed(seedSeq);
  }

  result_type operator() () {
    return engine();
  }

  bool operator == (const Impl& other) const {
    return engine == other.engine;
  }
};

// Engine's special member functions
Engine::Engine() : pImpl_(std::make_unique<Impl>()) {
#ifndef NDEBUG
  pImpl_->seed(272181374);
#else
  std::random_device randomDevice;

  pImpl_->seed(
    std::vector<int> {
      static_cast<int>(randomDevice()),
      static_cast<int>(randomDevice()),
      static_cast<int>(randomDevice()),
      static_cast<int>(randomDevice())
    }
  );
#endif
}
Engine::Engine(int seedArg) : Engine() {
  seed(seedArg);
}
Engine::Engine(Engine&& other) noexcept = default;
Engine& Engine::operator = (Engine&& other) noexcept = default;
Engine::Engine(const Engine& other) : pImpl_(std::make_unique<Impl>(*other.pImpl_)) {}
Engine& Engine::operator = (const Engine& other) {
  *pImpl_ = *other.pImpl_;
  return *this;
}
Engine::~Engine() = default;

void Engine::seed(int x) {
  pImpl_->seed(x);
}

void Engine::seed(const std::vector<int>& seeds) {
  pImpl_->seed(seeds);
}

Engine::result_type Engine::operator() () const {
  return pImpl_->operator()();
}

bool Engine::operator == (const Engine& other) const {
  return *pImpl_ == *other.pImpl_;
}

} // namespace Random
} // namespace Molassembler
} // namespace Scine
