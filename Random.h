#ifndef INCLUDE_TEMPLATE_MAGIC_RANDOM_H
#define INCLUDE_TEMPLATE_MAGIC_RANDOM_H

#include <random>
#include <vector>

/* TODO
 */

namespace TemplateMagic {

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
    randomEngine.seed(_seedSequence);
  }

public:
  mutable std::mt19937 randomEngine;

  Generator() {
    _initGenerator();
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
    std::vector<T> returnNumbers;
    std::uniform_real_distribution<T> uniformDistribution(lower, upper);

    for(unsigned i = 0; i < N; i++) {
      returnNumbers.emplace_back(uniformDistribution(randomEngine));
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
    std::vector<T> returnNumbers;
    std::uniform_int_distribution<T> uniformDistribution(lower, upper);

    for(unsigned i = 0; i < N; i++) {
      returnNumbers.emplace_back(uniformDistribution(randomEngine));
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
    std::uniform_real_distribution<T> uniformDistribution(lower, upper);
    return uniformDistribution(randomEngine);
  }

  template<typename T>
  std::enable_if_t<
    std::is_integral<T>::value,
    T
  > getSingle(
    const T& lower,
    const T& upper
  ) const {
    std::uniform_int_distribution<T> uniformDistribution(lower, upper);
    return uniformDistribution(randomEngine);
  }
};

static Generator random;

} // namespace TemplateMagic

#endif
