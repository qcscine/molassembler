// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include <boost/test/unit_test.hpp>

#include "temple/Optionals.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

struct Noisy {
  int f = -3;

  static unsigned nConstructed;

  Noisy() { ++nConstructed; }
  Noisy(const Noisy& /* other */) { ++nConstructed; }
  Noisy(Noisy&& /* other */) noexcept { ++nConstructed; }
  Noisy& operator = (const Noisy& other) = delete;
  Noisy& operator = (Noisy&& other) = delete;
  ~Noisy() = default;
};

unsigned Noisy::nConstructed = 0;

boost::optional<double> safe_root(double x) {
  if(x >= 0) {
    return {std::sqrt(x)};
  }

  return {};
}

boost::optional<double> safe_reciprocal(double x) {
  if(x != 0) {
    return 1.0 / x;
  }

  return {};
}

boost::optional<double> safe_binary(double x, double y) {
  double diff = y - x;
  if(diff != 0.0) {
    return 1 / diff;
  }

  return boost::none;
}

boost::optional<double> identity(double x) {
  return {x};
}

boost::optional<double> negate(double x) {
  return -x;
}

boost::optional<double> dWithCRefNoisy(double previous, const Noisy& x) {
  return previous + x.f;
}

boost::optional<int> withCRefNoisy(int previous, const Noisy& x) {
  return previous + x.f;
}

boost::optional<int> divFourBy(int previous) {
  if(previous != 0) {
    return static_cast<int>(
      4.0 / previous
    );
  }

  return boost::none;
}

boost::optional<int> safe_int_binary(const int& x, const int& y) {
  int diff = y - x;
  if(diff != 0) {
    return 1.0 / diff;
  }

  return boost::none;
}

template<typename T>
std::ostream& operator << (std::ostream& os, const boost::optional<T>& valueOptional) {
  if(valueOptional) {
    os << "Some " << valueOptional.value();
  } else {
    os << "None";
  }

  return os;
}

BOOST_AUTO_TEST_CASE(optionalEnhancementTests) {
  using namespace temple;

  const Noisy foo;

  auto first = safe_root(-4.3)
    | callIfNone(safe_binary, 2.3, 2.3)
    | callIfNone(identity, 4.9)
    | callIfNone(dWithCRefNoisy, -4, foo);

  std::cout << first << std::endl;

  auto second = safe_root(16.0)
    | callIfSome(negate, ANS)
    | callIfSome(safe_binary, 1, ANS);

  std::cout << second << std::endl;

  const int x = -4;

  auto third = boost::optional<int> {4}
    | callIfSome(withCRefNoisy, ANS, foo)
    | callIfSome(divFourBy, ANS)
    | callIfSome(safe_int_binary, ANS, x);

  std::cout << third << std::endl;

  auto typesTest = callIfSome(withCRefNoisy, 4, foo);
  using typesTestTuple = typename decltype(typesTest)::TupleType;
  static_assert(
    std::is_same<
      std::tuple_element_t<0, typesTestTuple>,
      int
    >::value,
    "No"
  );
  static_assert(
    std::is_same<
      std::tuple_element_t<1, typesTestTuple>,
      const Noisy&
    >::value,
    "No"
  );

  BOOST_CHECK(Noisy::nConstructed == 1);

  // Can there be different types?
  auto fourth = boost::optional<int> {0}
    | callIfSome(divFourBy, ANS)
    | callIfSome(safe_root, -4.3);

  std::cout << fourth << std::endl;

  /*static_assert(
    std::is_same<
      decltype(fourth),
      detail::Optional<double>
    >::value,
    "No"
  );*/
}
