/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Operator consistency checks
 *
 * Contains a suite of constexpr logical operator consistency checks. These can
 * aid in the diagnosis of custom operator weak ordering inconsistencies and
 * serve as useful library tests.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_LOGICAL_OPERATOR_TESTS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_LOGICAL_OPERATOR_TESTS_H

#include "Molassembler/Temple/constexpr/Math.h"

#include <iostream>

namespace Scine {
namespace Molassembler {
namespace Temple {

//! For any two types, check consistency of their logical operators
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

//! Limited variant of testLogicalOperators
template<typename T>
constexpr bool testOperatorSmaller(const T& a, const T& b) {
  return Math::XOR(
    a < b,
    b < a,
    !(a < b) && !(b < a) // a != b expressed with < only
  );
}

namespace dynamic {

//! Dynamic explainer of inconsistencies
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
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
