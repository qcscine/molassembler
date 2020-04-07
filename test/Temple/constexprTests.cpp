/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Temple/Adaptors/Zip.h"
#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/Random.h"
#include "molassembler/Temple/Stringify.h"

#include "molassembler/Temple/constexpr/Array.h"
#include "molassembler/Temple/constexpr/Bitmask.h"
#include "molassembler/Temple/constexpr/Bitset.h"
#include "molassembler/Temple/constexpr/Containers.h"
#include "molassembler/Temple/constexpr/DynamicArray.h"
#include "molassembler/Temple/constexpr/DynamicMap.h"
#include "molassembler/Temple/constexpr/DynamicSet.h"
#include "molassembler/Temple/constexpr/FloatingPointComparison.h"
#include "molassembler/Temple/constexpr/Jsf.h"
#include "molassembler/Temple/constexpr/LogicalOperatorTests.h"
#include "molassembler/Temple/constexpr/Math.h"
#include "molassembler/Temple/constexpr/TupleType.h"
#include "molassembler/Temple/constexpr/TupleTypePairs.h"

#include <iostream>
#include <iomanip>

#include <boost/test/results_collector.hpp>

using namespace Scine;
extern temple::Generator<> generator;

inline bool lastTestPassed() {
  using namespace boost::unit_test;
  test_case::id_t id = framework::current_test_case().p_id;
  test_results rez = results_collector.results(id);
  return rez.passed();
}

using namespace std::string_literals;

namespace ArrayTests {

constexpr auto testArr = temple::Array<unsigned, 3> {4, 3, 5};

template<typename T, size_t size>
constexpr temple::Array<T, size> modifyArray(
  const temple::Array<T, size>& array
) {
  auto arrayCopy = array;
  inPlaceSwap(arrayCopy, 0, 1);
  return arrayCopy;
}

constexpr auto modf = modifyArray(testArr);

static_assert(
  modf == temple::Array<unsigned, 3> {3, 4, 5},
  "Swap doesn't work as expected"
);

template<size_t size>
constexpr void testIteration(const temple::Array<unsigned, size>& array) {
  for(const auto& element : array) {
    std::cout << element << std::endl;
  }
}

static_assert(temple::Math::factorial(5) == 120, "Factorial is incorrect");
static_assert(temple::Math::factorial(0) == 1, "Factorial is incorrect");

static_assert(
  std::is_same<
    decltype(temple::makeArray(4, 3, 9)),
    temple::Array<int, 3>
  >::value,
  "makeArray does not work as expected"
);

} // namespace ArrayTests

BOOST_AUTO_TEST_CASE( mathApproxEqual ) {
  const unsigned numTests = 100;

  constexpr double accuracy = 1e-12;

  static_assert(
    accuracy >= std::numeric_limits<double>::epsilon(),
    "Testing accuracy must be greater than machine epsilon!"
  );

  // sqrt
  auto sqrt_failures = temple::copy_if(
    temple::random::getN<double>(0, 1e6, numTests, generator.engine),
    [&](const double randomPositiveNumber) -> bool {
      return !temple::floating::isCloseRelative(
        temple::Math::sqrt(randomPositiveNumber),
        std::sqrt(randomPositiveNumber),
        accuracy
      );
    }
  );

  BOOST_CHECK(sqrt_failures.empty());

  if(!sqrt_failures.empty()) {
    std::cout << "Square-root implementation is lacking! Failures: " << std::endl;

    for(const double value : sqrt_failures) {
      std::cout << "  x = " << std::setw(12) << value
        << ", sqrt = " << std::setw(12) << temple::Math::sqrt(value)
        << ", std::sqrt = " << std::setw(12) << std::sqrt(value)
        << ", |Δ| = " << std::setw(12) << std::fabs(
          temple::Math::sqrt(value) - std::sqrt(value)
        ) << std::endl;
    }

    std::cout << std::endl;
  }

  // asin
  const auto randomInverseTrigNumbers = temple::random::getN<double>(
    -1 + std::numeric_limits<double>::epsilon(),
    1 - std::numeric_limits<double>::epsilon(),
    numTests,
    generator.engine
  );

  auto asin_failures = temple::copy_if(
    randomInverseTrigNumbers,
    [&](const double randomInverseTrigNumber) -> bool {
      return !temple::floating::isCloseRelative(
        temple::Math::asin(randomInverseTrigNumber),
        std::asin(randomInverseTrigNumber),
        1e-8
      );
    }
  );

  BOOST_CHECK(asin_failures.empty());

  if(!asin_failures.empty()) {
    std::cout << "Inverse sine implementation is lacking! Failures: " << std::endl;

    for(const double value : asin_failures) {
      std::cout << "  x = " << std::setw(12) << value
        << ", asin = " << std::setw(12) << temple::Math::asin(value)
        << ", std::asin = " << std::setw(12) << std::asin(value)
        << ", |Δ| = " << std::setw(12) << std::fabs(
          temple::Math::asin(value) - std::asin(value)
        ) << std::endl;
    }

    std::cout << std::endl;
  }

  auto testPow = [&](const double number, const int& exponent) -> bool {
    const double test = temple::Math::pow(number, exponent);
    const double reference = std::pow(number, exponent);

    bool passes = temple::floating::isCloseRelative(
      test,
      reference,
      accuracy
    );

    if(!passes) {
      std::cout << "  x = " << std::setw(12) << number
        << ", exp = " << std::setw(4) << exponent
        << ", pow = " << std::setw(12) << test
        << ", std::pow = " << std::setw(12) << reference
        << ", |Δ| = " << std::setw(12) << std::fabs(test - reference) << ", max permissible diff: "
        << (
          accuracy * std::max(
            std::fabs(test),
            std::fabs(reference)
          )
        ) << std::endl;
    }

    return passes;
  };

  BOOST_CHECK(
    temple::all_of(
      temple::adaptors::zip(
        temple::random::getN<double>(-1e5, 1e5, numTests, generator.engine),
        temple::random::getN<int>(-40, 40, numTests, generator.engine)
      ),
      testPow
    )
  );

  auto testRecPow = [&](const double number, const unsigned exponent) -> bool {
    const double test = temple::Math::recPow(number, exponent);
    const double reference = std::pow(number, exponent);

    bool passes = temple::floating::isCloseRelative(
      test,
      reference,
      accuracy
    );

    if(!passes) {
      std::cout << "  x = " << std::setw(12) << number
        << ", exp = " << std::setw(4) << exponent
        << ", recPow = " << std::setw(12) << test
        << ", std::pow = " << std::setw(12) << reference
        << ", |Δ| = " << std::setw(12) << std::fabs(test - reference) << ", max permissible diff: "
        << (
          accuracy * std::max(
            std::fabs(test),
            std::fabs(reference)
          )
        ) << std::endl;
    }

    return passes;
  };

  BOOST_CHECK(
    temple::all_of(
      temple::adaptors::zip(
        temple::random::getN<double>(-1e5, 1e5, numTests, generator.engine),
        temple::random::getN<unsigned>(0, 40, numTests, generator.engine)
      ),
      testRecPow
    )
  );


  // ln
  const auto randomZ = temple::random::getN<double>(1e-10, 1e10, numTests, generator.engine);
  bool all_ln_pass = temple::all_of(
    randomZ,
    [&](const auto& z) -> bool {
      bool pass = temple::floating::isCloseRelative(
        temple::Math::ln(z),
        std::log(z),
        accuracy
      );

      if(!pass) {
        std::cout << "ln deviates for z = " << std::setw(12) << z
          << ", ln(z) = " << std::setw(12) << temple::Math::ln(z)
          << ", std::log(z) = " << std::setw(12) << std::log(z)
          << ", |Δ| = " << std::setw(12) << std::fabs(
            temple::Math::ln(z) - std::log(z)
          ) << std::endl;
      }

      return pass;
    }
  );

  BOOST_CHECK(all_ln_pass);

  BOOST_CHECK(
    temple::all_of(
      temple::random::getN<double>(-100, 100, numTests, generator.engine),
      [](const double x) -> bool {
        return(temple::Math::floor(x) <= x);
      }
    )
  );

  BOOST_CHECK(
    temple::all_of(
      temple::random::getN<double>(-100, 100, numTests, generator.engine),
      [](const double x) -> bool {
        return(temple::Math::ceil(x) >= x);
      }
    )
  );

  BOOST_CHECK(
    temple::all_of(
      temple::random::getN<double>(-M_PI / 2, M_PI / 2, numTests, generator.engine),
      [&](const double x) -> bool {
        return temple::floating::isCloseRelative(
          temple::Math::atan(x),
          std::atan(x),
          accuracy
        );
      }
    )
  );
}

BOOST_AUTO_TEST_CASE(arrayPermutation) {
  std::array<unsigned, 4> base {{0, 1, 2, 3}};
  std::array<unsigned, 4> STLComparison {{0, 1, 2, 3}};

  bool customHasNext, STLHasNext;

  do {
    customHasNext = temple::inPlaceNextPermutation(base);
    STLHasNext = std::next_permutation(STLComparison.begin(), STLComparison.end());

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In forward permutation, base is {" << temple::condense(base)
      << "} and and STL is {" << temple::condense(STLComparison)
      << "}"
    );
  } while(customHasNext && STLHasNext);

  BOOST_CHECK_MESSAGE(
   !customHasNext && !STLHasNext,
   "The two permutation algorithms don't terminate at the same time"
  );

  base = {{3, 2, 1, 0}};
  STLComparison = {{3, 2, 1, 0}};

  do {
    customHasNext = temple::inPlacePreviousPermutation(base);
    STLHasNext = std::prev_permutation(STLComparison.begin(), STLComparison.end());

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In backward permutation, base is {" << temple::condense(base)
      << "} and and STL is {" << temple::condense(STLComparison)
      << "}"
    );
  } while(customHasNext && STLHasNext);

  BOOST_CHECK_MESSAGE(
   !customHasNext && !STLHasNext,
   "The two permutation algorithms don't terminate at the same time"
  );

  // Variants of limited sections of the array

  base = {{0, 1, 2, 3}};
  STLComparison = {{0, 1, 2, 3}};

  do {
    customHasNext = temple::inPlaceNextPermutation(base, 1, 3);
    STLHasNext = std::next_permutation(STLComparison.begin() + 1, STLComparison.end() - 1);

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In limited forward permutation, base is {" << temple::condense(base)
      << "} and and STL is {" << temple::condense(STLComparison)
      << "}"
    );
  } while(customHasNext && STLHasNext);

  BOOST_CHECK_MESSAGE(
   !customHasNext && !STLHasNext,
   "The two permutation algorithms don't terminate at the same time in limited "
   "forward permutation"
  );

  base = {{3, 2, 1, 0}};
  STLComparison = {{3, 2, 1, 0}};

  do {
    customHasNext = temple::inPlacePreviousPermutation(base, 1, 3);
    STLHasNext = std::prev_permutation(STLComparison.begin() + 1, STLComparison.end() - 1);

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In limited backward permutation, base is {" << temple::condense(base)
      << "} and and STL is {" << temple::condense(STLComparison)
      << "}"
    );
  } while(customHasNext && STLHasNext);

  BOOST_CHECK_MESSAGE(
   !customHasNext && !STLHasNext,
   "The two permutation algorithms don't terminate at the same time in limited "
   "backward permutation"
  );
}

constexpr bool compileTimeDynTest() {
  temple::DynamicArray<unsigned, 10> nonConstArr {4, 3, 6};
  nonConstArr.push_back(9);
  return nonConstArr.size() == 4;
}

constexpr bool dynArrSpliceTest() {
  temple::DynamicArray<unsigned, 10> nonConstArr {4, 3, 6, 5, 1, 9};
  auto spliced = nonConstArr.splice(2);

  return (
    spliced == temple::DynamicArray<unsigned, 10> {6, 5, 1, 9}
    && nonConstArr == temple::DynamicArray<unsigned, 10> {4, 3}
  );
}

BOOST_AUTO_TEST_CASE(dynamicArrayTests) {
  constexpr temple::DynamicArray<unsigned, 10> arr {4, 3, 5};

  static_assert(
    arr.size() == 3,
    "Array size isn't initialized correctly from parameter pack ctor"
  );
  static_assert(
    compileTimeDynTest(),
    "non-const dynamic array push_back does not work as expected"
  );
  static_assert(
    dynArrSpliceTest(),
    "non-const dynamic array splice does not work as expected"
  );
  static_assert(
    arr.end() - arr.begin() == 3,
    "Subtracting begin/end iterators does not yield dynamic length"
  );

  static_assert(
    arr.begin() - arr.end() == -3,
    "Subtracting begin/end iterators does not yield dynamic length"
  );

  constexpr temple::Array<unsigned, 10> values {1, 2, 2, 3, 3, 3, 4, 4, 4, 4};

  constexpr auto grouped = groupByEquality(
    values,
    std::equal_to<>()
  );

  static_assert(
    grouped.size() == 4
    && grouped.at(0).size() == 1
    && grouped.at(1).size() == 2
    && grouped.at(2).size() == 3
    && grouped.at(3).size() == 4,
    "Grouping does not work as expected"
  );

  constexpr temple::DynamicArray<unsigned, 14> fromFixed {values};

  static_assert(fromFixed.size() == 10, "Construction from fixed doesn't work");
}

template<typename T, size_t size, class Comparator>
constexpr bool isSorted(const temple::DynamicSet<T, size, Comparator>& set) {
  Comparator comparator;

  auto left = set.begin();
  auto right = set.begin(); ++right;

  while(right != set.end()) {
    if(comparator(*right, *left)) {
      return false;
    }

    ++right;
    ++left;
  }

  return true;
}

BOOST_AUTO_TEST_CASE(dynamicSetTests) {
  temple::DynamicSet<unsigned, 10> set;

  BOOST_CHECK(set.size() == 0);
  BOOST_CHECK(
    std::distance(
      set.begin(),
      set.end()
    ) == 0u
  );

  for(const auto& item : {9u, 3u, 5u}) {
    set.insert(item);
  }

  BOOST_CHECK(set.size() == 3);
  BOOST_CHECK(
    std::distance(
      set.begin(),
      set.end()
    ) == 3u
  );
  BOOST_CHECK(set.contains(3) && set.contains(5) && set.contains(9));
  for(const auto& item : {2u, 4u, 8u, 10u}) {
    BOOST_CHECK_MESSAGE(
      !set.contains(item),
      "Set says it contains " << item << " when it shouldn't (set is {"
        << temple::condense(set)
        << "}."
    );
  }
  BOOST_CHECK(isSorted(set));

  temple::DynamicSet<unsigned, 10> setInitList {
    temple::DynamicArray<unsigned, 10> {
      4u, 9u, 13u
    }
  };

  BOOST_CHECK(setInitList.size() == 3);
  BOOST_CHECK(
    std::distance(
      setInitList.begin(),
      setInitList.end()
    ) == 3u
  );

  setInitList.insert(0u);

  BOOST_CHECK(setInitList.size() == 4);

  BOOST_CHECK_MESSAGE(
    setInitList.contains(4)
    && setInitList.contains(9)
    && setInitList.contains(13)
    && setInitList.contains(0)
    && !setInitList.contains(1)
    && !setInitList.contains(25),
    "setInitList {"
      << temple::condense(setInitList)
      << "} does not conform to expectations concerning contains:\n"
      << std::boolalpha
      << "contains 4, expect true:" << setInitList.contains(4)
      << "\ncontains 9, expect true:" << setInitList.contains(9)
      << "\ncontains 13, expect true:" << setInitList.contains(13)
      << "\ncontains 0, expect true:" << setInitList.contains(0)
      << "\ncontains 1, expect false:" << setInitList.contains(1)
      << "\ncontains 25, expect false:" << setInitList.contains(25)
  );
  // Set of arrays
  temple::Array<
    temple::Array<unsigned, 4>,
    5
  > sampleArrays {
    temple::Array<unsigned, 4> {1u, 2u, 3u, 4u},
    temple::Array<unsigned, 4> {1u, 2u, 4u, 3u},
    temple::Array<unsigned, 4> {1u, 4u, 3u, 2u},
    temple::Array<unsigned, 4> {1u, 4u, 2u, 3u},
    temple::Array<unsigned, 4> {2u, 1u, 3u, 4u}
  };

  temple::DynamicSet<
    temple::Array<unsigned, 4>,
    10
  > arraysSet;

  arraysSet.insert(sampleArrays.at(0));
  arraysSet.insert(sampleArrays.at(2));
  arraysSet.insert(sampleArrays.at(3));

  BOOST_CHECK(arraysSet.contains(sampleArrays.at(0)));
  BOOST_CHECK(arraysSet.contains(sampleArrays.at(2)));
  BOOST_CHECK(arraysSet.contains(sampleArrays.at(3)));
  BOOST_CHECK(!arraysSet.contains(sampleArrays.at(1)));
  BOOST_CHECK(!arraysSet.contains(sampleArrays.at(4)));

  BOOST_CHECK(arraysSet.size() == 3);
  BOOST_CHECK(
    std::distance(
      arraysSet.begin(),
      arraysSet.end()
    ) == 3
  );
}

template<typename T, size_t size>
bool validate(const temple::DynamicSet<T, size>& set) {
  // Is the set ordered?
  auto leftIter = set.begin();
  auto rightIter = leftIter; ++rightIter;
  auto bound = set.end(); --bound;

  if(leftIter == set.end() || rightIter == set.end()) {
    return true;
  }

  while(rightIter != bound) {
    if(*leftIter >= *rightIter) {
      std::cout << "*left >= *right -> " << *leftIter << " >= " << *rightIter << std::endl;
      return false;
    }

    ++leftIter;
    ++rightIter;
  }

  // Is the reported size equal to a begin-end through-iteration?
  return (
    set.size() == static_cast<unsigned>(
      std::distance(
        set.begin(),
        set.end()
      )
    )
  );
}

BOOST_AUTO_TEST_CASE(arrayOperators) {
  temple::Array<unsigned, 4> a {4, 2, 3, 1};
  temple::Array<unsigned, 4> b {4, 3, 2, 1};

  BOOST_CHECK(temple::testLogicalOperators(a, b));
  BOOST_CHECK(temple::testLogicalOperators(a, a));

  temple::dynamic::explainLogicalOperatorFailures(a, b);
}

BOOST_AUTO_TEST_CASE(dynamicSetFuzzing) {
  for(unsigned N = 0; N < 100; ++N) {
    temple::DynamicSet<unsigned, 100> subject;

    std::vector<unsigned> numbers;
    numbers.resize(50);
    std::iota(
      numbers.begin(),
      numbers.end(),
      0
    );

    temple::random::shuffle(numbers, generator.engine);

    for(const auto& number : numbers) {
      subject.insert(number);

      bool isValid = validate(subject);
      if(!isValid) {
        std::cout << "After inserting " << number
          << ", set is left in invalid state. Set: {"
          << temple::condense(subject) << "}\ninsert sequence {"
          << temple::condense(numbers) << "}."
          << std::endl;
      }

      BOOST_REQUIRE(isValid);

      BOOST_CHECK(subject.contains(number));
    }
  }
}

namespace tuple_type_tests {

struct Apple {
  static constexpr unsigned number = 4;
};

struct Banana {
  static constexpr unsigned number = 3;
};

struct Cherry {
  static constexpr unsigned number = 50;
};

using Fruit = std::tuple<Apple, Banana, Cherry>;

template<typename ... FruitClasses>
struct sumNumbersValue {
  static constexpr unsigned initializer() {
    const std::array<unsigned, sizeof...(FruitClasses)> sizes = {{
      FruitClasses::number...
    }};

    unsigned sum = 0;

    for(unsigned i = 0; i < sizes.size(); ++i) {
      sum += sizes[i];
    }

    return sum;
  };

  static constexpr unsigned value = initializer();
};

template<typename ... FruitClasses>
struct sumNumbersFunctor {
  static constexpr unsigned value() {
    const std::array<unsigned, sizeof...(FruitClasses)> sizes = {{
      FruitClasses::number...
    }};

    unsigned sum = 0;

    for(unsigned i = 0; i < sizes.size(); ++i) {
      sum += sizes[i];
    }

    return sum;
  };
};

template<typename FruitClass>
struct getNumberFunctor {
  static constexpr unsigned value() {
    return FruitClass::number;
  }
};

template<typename FruitClass>
struct getNumberValue {
  static constexpr unsigned initialize() {
    return FruitClass::number;
  }

  static constexpr unsigned value = initialize();
};

template<typename FruitClass>
constexpr unsigned getNumberFunction() {
  return FruitClass::number;
}

template<typename A, typename B>
struct pairSumFunctor {
  static constexpr unsigned value() {
    return A::number + B::number;
  }
};

template<typename A, typename B>
struct pairSumValue {
  static constexpr unsigned initialize() {
    return A::number + B::number;
  }

  static constexpr unsigned value = initialize();
};

static_assert(
  temple::tuples::unpackToFunction<Fruit, sumNumbersValue>() == 57,
  "Unpacking fruit tuple to valueValue does not yield expected result"
);

static_assert(
  temple::tuples::unpackToFunction<Fruit, sumNumbersFunctor>() == 57,
  "Unpacking fruit tuple to valueFunctor does not yield expected result"
);

static_assert(
  temple::arraysEqual(
    temple::tuples::map<Fruit, getNumberValue>(),
    temple::Array<unsigned, 3> {{Apple::number, Banana::number, Cherry::number}}
  ),
  "Mapping with getNumberValue does not yield expected result!"
);

static_assert(
  temple::arraysEqual(
    temple::tuples::map<Fruit, getNumberFunctor>(),
    temple::Array<unsigned, 3> {{Apple::number, Banana::number, Cherry::number}}
  ),
  "Mapping with getNumberFunctor does not yield expected result!"
);

static_assert(
  temple::arraysEqual(
    temple::tuples::mapAllPairs<Fruit, pairSumValue>(),
    temple::Array<unsigned, 3> {{
      Apple::number + Banana::number,
      Apple::number + Cherry::number,
      Banana::number + Cherry::number
    }}
  ),
  "Mapping pairs with pairSumValue does not yield expected result!"
);

static_assert(
  temple::arraysEqual(
    temple::tuples::mapAllPairs<Fruit, pairSumFunctor>(),
    temple::Array<unsigned, 3> {{
      Apple::number + Banana::number,
      Apple::number + Cherry::number,
      Banana::number + Cherry::number
    }}
  ),
  "Mapping pairs with pairSumFunctor does not yield expected result!"
);

using countTestType = std::tuple<unsigned, float, double, unsigned, size_t>;

static_assert(
  temple::tuples::countType<countTestType, unsigned>() == 2,
  "Counting unsigned in countTestType does not return two!"
);

} // namespace tuple_type_tests

namespace FloatingPointComparisonTests {

template<typename T>
constexpr bool testAbsoluteComparison(const T& a, const T& b, const T& tolerance) {
  temple::floating::ExpandedAbsoluteEqualityComparator<T> comparator {
    tolerance
  };

  return (
    temple::Math::XOR(
      (
        comparator.isLessThan(a, b)
        && comparator.isMoreThan(b, a)
        && comparator.isUnequal(a, b)
      ),
      (
        comparator.isLessThan(b, a)
        && comparator.isMoreThan(a, b)
        && comparator.isUnequal(a, b)
      ),
      (
        !comparator.isLessThan(a, b)
        && !comparator.isMoreThan(a, b)
        && comparator.isEqual(a, b)
      )
    ) && temple::Math::XOR(
      comparator.isEqual(a, b),
      comparator.isUnequal(a, b)
    )
  );
}

template<typename T>
constexpr bool testRelativeComparison(const T& a, const T& b, const T& tolerance) {
  temple::floating::ExpandedRelativeEqualityComparator<T> comparator {
    tolerance
  };

  return (
    temple::Math::XOR(
      (
        comparator.isLessThan(a, b)
        && comparator.isMoreThan(b, a)
        && comparator.isUnequal(a, b)
      ),
      (
        comparator.isLessThan(b, a)
        && comparator.isMoreThan(a, b)
        && comparator.isUnequal(a, b)
      ),
      (
        !comparator.isLessThan(a, b)
        && !comparator.isMoreThan(a, b)
        && comparator.isEqual(a, b)
      )
    ) && temple::Math::XOR(
      comparator.isEqual(a, b),
      comparator.isUnequal(a, b)
    )
  );
}

static_assert(
  testAbsoluteComparison(4.3, 3.9, 1e-4)
  && testAbsoluteComparison(4.3, 3.9, 1.0)
  && testAbsoluteComparison(4.4, 4.4, 1e-10),
  "absolute comparison has inconsistent operators!"
);

static_assert(
  testRelativeComparison(4.3, 3.9, 1e-4)
  && testRelativeComparison(4.3, 3.9, 1.0)
  && testRelativeComparison(4.4, 4.4, 1e-10),
  "relative comparison has inconsistent operators!"
);

} // namespace FloatingPointComparisonTests

namespace ConcatenationTests {

constexpr temple::Array<unsigned, 4> f {4, 2, 9, 3};
constexpr temple::Array<unsigned, 4> g {11, 22, 33, 44};
constexpr temple::Array<unsigned, 4> h {234, 292, 912, 304};
constexpr temple::Array<unsigned, 8> fg {
  4, 2, 9, 3,
  11, 22, 33, 44
};
constexpr temple::Array<unsigned, 12> fgh {
  4, 2, 9, 3,
  11, 22, 33, 44,
  234, 292, 912, 304
};

static_assert(
  temple::arraysEqual(
    temple::arrayConcatenate(f, g),
    fg
  ),
  "Pairwise concatenation does not preserve sequence!"
);

static_assert(
  temple::arraysEqual(
    temple::arrayConcatenate(f, g, h),
    fgh
  ),
  "Variadic concatenation does not work as expected"
);

} // namespace ConcatenationTests

namespace DynamicMapTests {

constexpr temple::DynamicMap<unsigned, int, 20> generateMap() {
  temple::DynamicMap<unsigned, int, 20> myMap;

  myMap.insert(4, -2);
  myMap.insert(1, 4);
  myMap.insert(3, 9);

  return myMap;
}

constexpr auto a = generateMap();

static_assert(a.at(4u) == -2, "Map does not find element with key 4");
static_assert(a.at(1u) == 4, "Map does not find element with key 1");
static_assert(a.at(3u) == 9, "Map does not find element with key 3");

} // namespace DynamicMapTests

namespace UpperTriangularMatrixTests {

// Can default-construct
constexpr auto defaultMatr = temple::UpperTriangularMatrix<bool, 15> {};
static_assert(decltype(defaultMatr)::N == 6u, "Size isn't right");

constexpr auto matr = temple::makeUpperTriangularMatrix(
  std::array<unsigned, 6> {{1, 2, 3, 4, 5, 6}}
);

static_assert(decltype(matr)::N == 4u, "Size isn't right");

/*constexpr auto failing = temple::makeUpperTriangularMatrix(
  std::array<unsigned, 5> {{1, 2, 3, 4, 5}}
);*/

constexpr auto fromArray = temple::makeUpperTriangularMatrix(
  temple::Array<unsigned, 6> {{1, 2, 3, 4, 5, 6}}
);
static_assert(decltype(fromArray)::N == 4u, "Size isn't right");

} // namespace UpperTriangularMatrixTests

enum class ScopedEnum : unsigned {A, B, C};
enum UnscopedEnum : unsigned {D, E};

BOOST_AUTO_TEST_CASE(bitmaskAll) {
  using namespace temple;

  constexpr auto a = make_bitmask(ScopedEnum::A) | ScopedEnum::C;
  static_assert(a & ScopedEnum::A, "A must be contained in the bitmask a");
  BOOST_CHECK(a[ScopedEnum::A] && a[ScopedEnum::C] && !a[ScopedEnum::B]);

  constexpr auto b = make_bitmask(UnscopedEnum::E);
  BOOST_CHECK(b[UnscopedEnum::E] && !b[UnscopedEnum::D]);
}

constexpr temple::Bitset<300> make_a_bitset() {
  temple::Bitset<300> a;

  a.set(0);
  a.set(299);
  a.set(63);
  a.set(128);

  return a;
}

BOOST_AUTO_TEST_CASE(bitsetTests) {
  using namespace temple;

  constexpr auto a = make_a_bitset();

  static_assert(a.test(0), "Zero hasn't been set");

  BOOST_CHECK(a.test(0) && a.test(299) && a.test(63) && a.test(128));
  BOOST_CHECK(!a.test(1) && !a.test(298) && !a.test(64) && !a.test(127));

  auto b = make_a_bitset();
  b.unset(63);
  BOOST_CHECK(!b.test(63));
}

BOOST_AUTO_TEST_CASE(permutationIndexTests) {
  using namespace temple;

  auto a = iota<Array, unsigned, 6>();

  // expect monotonically increasing permutation index on inplaceNextPermutation
  auto index = permutationIndex(a);

  BOOST_CHECK(index == 0);

  bool pass = true;
  auto previousIndex = index;
  while(inPlaceNextPermutation(a)) {
    index = permutationIndex(a);
    if(index != previousIndex + 1) {
      pass = false;
      break;
    }
    previousIndex = index;
  }

  BOOST_CHECK_MESSAGE(
    pass,
    "inPlaceNextPermutation must yield monotonically increasing permutation "
    "index values"
  );
}
