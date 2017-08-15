#ifndef INCLUDE_CONSTEXPR_MAGIC_LOGICAL_OPERATOR_TESTS_H
#define INCLUDE_CONSTEXPR_MAGIC_LOGICAL_OPERATOR_TESTS_H

#include "Math.h"
#include <iostream>

namespace ConstexprMagic {

// For any two types, check consistency of their logical operators
template<typename T>
constexpr bool testLogicalOperators(const T& a, const T& b) {
  return (
    Math::XOR( // only one of the following three cases may be true at any time
      a < b && b > a && a != b,
      b < a && a > b && a != b,
      !(a < b) && !(a > b) && a == b
    ) && Math::XOR( // ensure == and != are correct
      a == b,
      a != b
    )
  );
}

template<typename T>
constexpr bool testOperatorSmaller(const T& a, const T& b) {
  return Math::XOR(
    a < b,
    b < a, 
    !(a < b) && !(b < a) // a != b expressed with < only
  );
}

namespace dynamic {

template<typename T>
void explainLogicalOperatorFailures(const T& a, const T& b) {
  if(
    !Math::XOR(
      a < b && b > a && a != b,
      b < a && a > b && a != b,
      !(a < b) && !(a > b) && a == b
    )
  ) {
    std::cout << "operator < is inconsistent:\n" << std::boolalpha
      << " a < b && b > a && a != b -> " 
      << (a < b) << " && " << (b > a) << " && " << (a != b) << " -> "
      << (a < b && b > a && a != b) << "\n"
      << " b < a && a > b && a != b -> "
      << (b < a) << " && " << (a > b) << " && " << (a != b) << " -> "
      << (b < a && a > b && a != b) << "\n"
      << " !(a < b) && !(a > b) && a == b -> "
      << !(a < b) << " && " << !(a > b) << " && " << (a == b) << " -> "
      << (!(a < b) && !(a > b) && a == b) << "\n";
  }

  if(
    !Math::XOR(
      a == b,
      a != b
    )
  ) {
    std::cout << "operator == is inconsistent:\n" << std::boolalpha
      << " a == b -> " << (a == b) << "\n"
      << " a != b -> " << (a != b) << "\n";
  }
}

} // namespace dynamic

} // namespace ConstexprMagic

#endif
