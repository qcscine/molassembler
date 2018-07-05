#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_RANDOM_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_RANDOM_H

#include "constexpr/JSF.h"
#include "Traits.h"

#include <random>

/*! @file
 *
 * Provides both an initialized randomness engine as well as random uniform
 * integers or floating-point numbers.
 */

namespace temple {

class Generator {
public:
  mutable jsf::JSF32 engine;

private:
  static std::array<uint32_t, 4> getRandomSeeds() {
    std::random_device randomDevice;

    return {
      randomDevice(),
      randomDevice(),
      randomDevice(),
      randomDevice()
    };
  }

  void _initializeEngine () {
#ifndef NDEBUG
    engine.seed(272181374);
#else
    engine.seed(getRandomSeeds());
#endif
  }

public:
  explicit Generator() {
    _initializeEngine();
  }

  template<typename Container>
  std::enable_if_t<
    std::is_signed<
      traits::getValueType<Container>
    >::value,
    void
  > seed(const Container& seeds) {
    std::seed_seq seedSeq {
      std::begin(seeds),
      std::end(seeds)
    };

    engine.seed(seedSeq);
  }

  void seed(int x) {
    engine.seed(x);
  }

  template<typename T>
  std::enable_if_t<
    std::is_floating_point<T>::value,
    std::vector<T>
  > getN(
    const T& lower,
    const T& upper,
    const unsigned N
  ) const {
    assert(lower <= upper);
    std::vector<T> returnNumbers;
    returnNumbers.reserve(N);

    std::uniform_real_distribution<T> uniformDistribution(lower, upper);

    for(unsigned i = 0; i < N; i++) {
      returnNumbers.emplace_back(uniformDistribution(engine));
    }

    return returnNumbers;
  }

  template<typename T>
  std::enable_if_t<
    std::is_integral<T>::value,
    std::vector<T>
  > getN(
    const T& lower,
    const T& upper,
    const unsigned N
  ) const {
    assert(lower <= upper);
    std::vector<T> returnNumbers;
    returnNumbers.reserve(N);

    std::uniform_int_distribution<T> uniformDistribution(lower, upper);

    for(unsigned i = 0; i < N; i++) {
      returnNumbers.emplace_back(uniformDistribution(engine));
    }

    return returnNumbers;
  }

  template<typename T>
  std::enable_if_t<
    std::is_floating_point<T>::value,
    T
  > getSingle(
    const T& lower,
    const T& upper
  ) const {
    assert(lower <= upper);
    std::uniform_real_distribution<T> uniformDistribution(lower, upper);
    return uniformDistribution(engine);
  }

  template<typename T>
  std::enable_if_t<
    std::is_integral<T>::value,
    T
  > getSingle(
    const T& lower,
    const T& upper
  ) const {
    assert(lower <= upper);
    std::uniform_int_distribution<T> uniformDistribution(lower, upper);
    return uniformDistribution(engine);
  }

  template<typename T>
  std::enable_if_t<
    std::is_same<T, bool>::value,
    bool
  > getSingle() const {
    std::uniform_int_distribution<unsigned> uniformDistribution(0, 1);
    return static_cast<bool>(
      uniformDistribution(engine)
    );
  }

  template<typename Container>
  unsigned pickDiscrete(const Container& weights) {
    std::discrete_distribution<unsigned> distribution {
      std::begin(weights),
      std::end(weights)
    };

    return distribution(engine);
  }

  template<typename Container>
  void shuffle(Container& container) {
    std::shuffle(
      std::begin(container),
      std::end(container),
      engine
    );
  }
};

} // namespace temple

#endif
