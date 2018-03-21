#ifndef INCLUDE_LOGICAL_OPERATOR_TESTS_H
#define INCLUDE_LOGICAL_OPERATOR_TESTS_H

namespace OperatorTests {

// n-argument boolean XOR with the definition: only one argument may be true
constexpr unsigned TMPSum() {
  return 0;
}

template<typename T1, typename... T>
constexpr unsigned TMPSum(T1 a, T ... pack) {
  return a + TMPSum(pack ...);
}

template<typename ... Bools>
constexpr bool XOR(Bools ... bools) {
  return TMPSum(bools ...) == 1;
}

// For any two types, check consistency of their logical operators
template<typename T>
bool testLogicalOperators(const T& a, const T& b) {
  return (
    XOR( // only one of the following three cases may be true at any time
      a < b && b > a && a != b,
      b < a && a > b && a != b,
      !(a < b) && !(a > b) && a == b
    ) && XOR( // ensure == and != are correct
      a == b,
      a != b
    )
  );
}

template<typename T>
bool testOperatorSmaller(const T& a, const T& b) {
  return XOR(
    a < b,
    b < a, 
    !(a < b) && !(b < a) // a != b expressed with < only
  );
}

} // namespace OperatorTests

#endif
