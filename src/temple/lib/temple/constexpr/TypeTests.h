// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_TYPE_TESTS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_TYPE_TESTS_H

#include "temple/constexpr/Math.h"

#include <iostream>

/*! @file
 *
 * @brief Logical operator test helpers
 *
 * Contains a suite of constexpr logical operator consistency checks. These can
 * aid in the diagnosis of custom operator weak ordering inconsistencies and
 * serve as useful library tests.
 */

namespace temple {

namespace TypeTests {

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

template<class Container>
constexpr bool testConstIterators(const Container& container) {
  if(!( // Very basic tests
    container.begin() == container.begin()
    && container.end() == container.end()
    && container.begin() != container.end()
  )) {
    return false;
  }

  return true;
}

template<class Container>
constexpr bool testIterators(Container container) {
  /* TODO
   * - different tests depending on the type of iterator
   *   -> forward_iterator / bidirectional_iterator / random_access_iterator
   * - Maybe with std::advance
   * - see tests.cpp for the array types, I think there are some RAIter tests
   *   there
   */

  if(!testConstIterators(container)) {
    return false;
  }

  if(!(
    container.begin() == container.begin()
    && container.end() == container.end()
    && container.begin() != container.end()
  )) {
    return false;
  }

  return true;
}

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

} // namespace TypeTests

} // namespace temple

#endif
