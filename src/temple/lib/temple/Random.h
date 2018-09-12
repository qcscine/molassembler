#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_RANDOM_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_RANDOM_H

#include "temple/constexpr/JSF.h"
#include "temple/Traits.h"

#include <algorithm>

/*! @file
 *
 * @brief Randomness helper functions
 *
 * Provides both an initialized randomness engine as well as random uniform
 * integers or floating-point numbers.
 */

namespace temple {

class Generator {
public:
//!@name State
//!@{
  //! Underyling PRNG engine
  mutable jsf::JSF32 engine;
//!@}

private:
  //! Generate random seeds from a hardware random device (if STL implemented)
  static std::array<uint32_t, 4> getRandomSeeds() {
    std::random_device randomDevice;

    return {
      randomDevice(),
      randomDevice(),
      randomDevice(),
      randomDevice()
    };
  }

  //! Initialize the underlying engine
  void _initializeEngine () {
#ifndef NDEBUG
    engine.seed(272181374);
#else
    engine.seed(getRandomSeeds());
#endif
  }

public:
  //! Default constructor
  explicit Generator() {
    _initializeEngine();
  }

//!@name Modification
//!@{
  //! Re-seed the underlying PRNG
  void seed(int x) {
    engine.seed(x);
  }

  //! Re-seed the underlying PRNG
  void seed(const std::vector<int>& signedSeeds) {
    std::seed_seq seedSeq(
      std::begin(signedSeeds),
      std::end(signedSeeds)
    );

    engine.seed(seedSeq);
  }
//!@}

//!@name Generators
//!@{
  //! Generate N floating point values in the range [lower, upper)
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

  //! Generate N integral values in the range [lower, upper]
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

  //! Generate a single floating point value in the range [lower, upper)
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

  //! Generate a single integral value in the range [lower, upper]
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

  //! Flip a coin
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

  //! Pick a value using provided discrete distribution weights
  template<class Container>
  unsigned pickDiscrete(const Container& weights) {
    std::discrete_distribution<unsigned> distribution {
      std::begin(weights),
      std::end(weights)
    };

    return distribution(engine);
  }
//!@}

//!@name Auxiliary
//!@{
  //! Use underlying PRNG to shuffle a Container in-place
  template<class Container>
  void shuffle(Container& container) {
    std::shuffle(
      std::begin(container),
      std::end(container),
      engine
    );
  }
//!@}
};

} // namespace temple

#endif
