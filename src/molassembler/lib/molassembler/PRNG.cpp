#include "PRNG.h"

#include "temple/constexpr/JSF.h"
#include <random>

namespace molassembler {

namespace random {

struct Engine::Impl {
  temple::jsf::JSF32 engine;

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
};

// Engine's special member functions
Engine::Engine() : _pImpl(std::make_unique<Impl>()) {
#ifndef NDEBUG
  _pImpl->seed(272181374);
#else
  std::random_device randomDevice;

  return {
    randomDevice(),
    randomDevice(),
    randomDevice(),
    randomDevice()
  };

  _pImpl->seed(getRandomSeeds());
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

// Engine's member functions
Engine::result_type Engine::min() {
  return 0;
}

Engine::result_type Engine::max() {
  return ~result_type(0);
}

void Engine::seed(int x) {
  _pImpl->seed(x);
}

void Engine::seed(const std::vector<int>& seeds) {
  _pImpl->seed(seeds);
}

Engine::result_type Engine::operator() () const {
  return _pImpl->operator()();
}

} // namespace random

} // namespace molassembler
