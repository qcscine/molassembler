/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "Prng.h"

#include "temple/constexpr/Jsf.h"
#include <random>

namespace Scine {
namespace molassembler {
namespace random {

struct Engine::Impl {
  temple::JSF32 engine;

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
Engine::Engine() : _pImpl(std::make_unique<Impl>()) {
#ifndef NDEBUG
  _pImpl->seed(272181374);
#else
  std::random_device randomDevice;

  _pImpl->seed(
    std::vector<int> {
      static_cast<int>(randomDevice()),
      static_cast<int>(randomDevice()),
      static_cast<int>(randomDevice()),
      static_cast<int>(randomDevice())
    }
  );
#endif
}
Engine::Engine(Engine&& other) noexcept = default;
Engine& Engine::operator = (Engine&& other) noexcept = default;
Engine::Engine(const Engine& other) : _pImpl(std::make_unique<Impl>(*other._pImpl)) {}
Engine& Engine::operator = (const Engine& other) {
  *_pImpl = *other._pImpl;
  return *this;
}
Engine::~Engine() = default;

void Engine::seed(int x) {
  _pImpl->seed(x);
}

void Engine::seed(const std::vector<int>& seeds) {
  _pImpl->seed(seeds);
}

Engine::result_type Engine::operator() () const {
  return _pImpl->operator()();
}

bool Engine::operator == (const Engine& other) const {
  return *_pImpl == *other._pImpl;
}

} // namespace random
} // namespace molassembler
} // namespace Scine
