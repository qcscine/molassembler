#ifndef INCLUDE_TEMPLATE_MAGIC_RANDOM_H
#define INCLUDE_TEMPLATE_MAGIC_RANDOM_H

#include "Traits.h"

#include <random>

/*! @file
 *
 * Provides both an initialized randomness engine as well as random uniform
 * integers or floating-point numbers.
 */

namespace temple {

class Generator {
private:
  std::vector<unsigned> _seeds;

  void _initGenerator() {

#ifdef NDEBUG
    std::random_device randomDevice;
    for(unsigned n = 0; n < 5; n++) _seeds.emplace_back(randomDevice());
#else
    _seeds.emplace_back(2721813754);
#endif

    std::seed_seq _seedSequence(_seeds.begin(), _seeds.end());
    engine.seed(_seedSequence);
  }

public:
  mutable std::mt19937 engine;

  template<typename Container>
  std::enable_if_t<
    std::is_same<
      traits::getValueType<Container>,
      unsigned
    >::value,
    void
  > seed(const Container& container) {
    _seeds.clear();

    std::copy(
      std::begin(container),
      std::end(container),
      std::back_inserter(_seeds)
    );

    std::seed_seq seedSequence {
      std::begin(_seeds),
      std::end(_seeds)
    };

    engine.seed(seedSequence);
  }


  Generator() {
    _initGenerator();
  }

  template<typename Container>
  explicit Generator(const Container& seeds) {
    seed(seeds);
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

  const std::vector<unsigned>& getSeeds() const {
    return _seeds;
  }
};

} // namespace temple

#endif
