#define BOOST_TEST_MODULE ConstexprMagicTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "template_magic/Containers.h"
#include "template_magic/Random.h"

#include "Math.h"
#include "Containers.h"
#include "Array.h"
#include "Set.h"
#include "DynamicArray.h"
#include "DynamicSet.h"
#include "TupleType.h"
#include "LogicalOperatorTests.h"
#include "FloatingPointComparison.h"

#include <iostream>
#include <iomanip>

using namespace std::string_literals;

namespace ArrayTests {

constexpr auto testArr = ConstexprMagic::Array<unsigned, 3> {4, 3, 5};

template<typename T, size_t size>
constexpr ConstexprMagic::Array<T, size> modifyArray(
  const ConstexprMagic::Array<T, size>& array
) {
  auto arrayCopy = array;
  inPlaceSwap(arrayCopy, 0, 1);
  return arrayCopy;
}
  
constexpr auto modf = modifyArray(testArr);

static_assert(
  modf == ConstexprMagic::Array<unsigned, 3> {3, 4, 5},
  "Swap doesn't work as expected"
);

static_assert(
  ConstexprMagic::arrayPop(testArr) == ConstexprMagic::Array<unsigned, 2> {4, 3},
  "Pop doesn't work"
);

static_assert(
  ConstexprMagic::arrayPush(testArr, 9u) == ConstexprMagic::Array<unsigned, 4> {4, 3, 5, 9},
  "Push doesn't work"
);

static_assert(
  ConstexprMagic::arrayPush(testArr, 9u) == ConstexprMagic::Array<unsigned, 4> {4, 3, 5, 9},
  "arrayPush doesn't work on ConstexprMagic::Array"
);

constexpr auto stdTestArr = std::array<unsigned, 3> {4, 3, 5};

static_assert(
  ConstexprMagic::arraysEqual(ConstexprMagic::arrayPush(stdTestArr, 9u), std::array<unsigned, 4> {4, 3, 5, 9}),
  "arrayPush doesn't work on std::array"
);

template<size_t size>
constexpr void testIteration(const ConstexprMagic::Array<unsigned, size>& array) {
  for(const auto& element : array) {
    std::cout << element << std::endl;
  }
}

constexpr auto sortedArr = ConstexprMagic::Array<unsigned, 4> {4, 6, 9, 11};
constexpr auto oneMore = ConstexprMagic::insertIntoSorted(sortedArr, 5u);

static_assert(
  oneMore == ConstexprMagic::Array<unsigned, 5> {4, 5, 6, 9, 11},
  "InsertIntoSorted does not work as expected."
);

static_assert(ConstexprMagic::Math::factorial(5) == 120, "Factorial is incorrect");

// C++17 with std::array
/*
constexpr auto stdSortedArr = std::array<unsigned, 4> {{4, 6, 9, 11}};
constexpr auto stdOneMore = ConstexprMagic::insertIntoSorted(stdSortedArr, 5u);

static_assert(
  ConstexprMagic::arraysEqual(stdOneMore, std::array<unsigned, 5> {4, 5, 6, 9, 11}),
  "InsertIntoSorted does not work as expected with std::array"
);
// std::array::operator == (const std::array& other) isn't constexpr in C++17
*/

constexpr auto testSet = ConstexprMagic::makeSetFromSortedArray(oneMore);
constexpr auto withElement = testSet.insert(5);

} // namespace ArrayTests

BOOST_AUTO_TEST_CASE( mathApproxEqual ) {
  const unsigned numTests = 100;

  constexpr double accuracy = 1e-12;

  static_assert(
    accuracy >= std::numeric_limits<double>::epsilon(),
    "Testing accuracy must be greater than machine epsilon!"
  );

  // sqrt
  const auto randomPositiveNumbers = TemplateMagic::random.getN<double>(0, 1e6, numTests);
  auto sqrt_passes = TemplateMagic::map(
    randomPositiveNumbers,
    [](const double& randomPositiveNumber) -> bool {
      return ConstexprMagic::floating::isCloseRelative(
        ConstexprMagic::Math::sqrt(randomPositiveNumber),
        std::sqrt(randomPositiveNumber),
        accuracy
      );
    }
  );

  bool all_sqrt_pass = TemplateMagic::all_of(sqrt_passes);

  BOOST_CHECK(all_sqrt_pass);

  if(!all_sqrt_pass) {
    std::cout << "Square-root implementation is lacking! Failures: " << std::endl;

    for(unsigned i = 0; i < sqrt_passes.size(); i++) {
      if(!sqrt_passes[i]) {
        std::cout << "  x = " << std::setw(12) << randomPositiveNumbers[i]
          << ", sqrt = " << std::setw(12) << ConstexprMagic::Math::sqrt(randomPositiveNumbers[i])
          << ", std::sqrt = " << std::setw(12) << std::sqrt(randomPositiveNumbers[i])
          << ", |Δ| = " << std::setw(12) << std::fabs(
            ConstexprMagic::Math::sqrt(randomPositiveNumbers[i])
            - std::sqrt(randomPositiveNumbers[i])
          ) << std::endl;
      }
    }

    std::cout << std::endl;
  }

  // asin
  const auto randomInverseTrigNumbers = TemplateMagic::random.getN<double>(
    -1 + std::numeric_limits<double>::epsilon(),
    1 - std::numeric_limits<double>::epsilon(),
    numTests
  );

  auto asin_passes = TemplateMagic::map(
    randomInverseTrigNumbers,
    [](const double& randomInverseTrigNumber) -> bool {
      return ConstexprMagic::floating::isCloseRelative(
        ConstexprMagic::Math::asin(randomInverseTrigNumber),
        std::asin(randomInverseTrigNumber),
        accuracy
      );
    }
  );

  bool all_asin_pass = TemplateMagic::all_of(asin_passes);

  BOOST_CHECK(all_asin_pass);

  if(!all_asin_pass) {
    std::cout << "Inverse sine implementation is lacking! Failures: " << std::endl;

    for(unsigned i = 0; i < asin_passes.size(); i++) {
      if(!asin_passes[i]) {
        std::cout << "  x = " << std::setw(12) << randomInverseTrigNumbers[i]
          << ", asin = " << std::setw(12) << ConstexprMagic::Math::asin(randomInverseTrigNumbers[i])
          << ", std::asin = " << std::setw(12) << std::asin(randomInverseTrigNumbers[i])
          << ", |Δ| = " << std::setw(12) << std::fabs(
            ConstexprMagic::Math::asin(randomInverseTrigNumbers[i])
            - std::asin(randomInverseTrigNumbers[i])
          ) << std::endl;
      }
    }

    std::cout << std::endl;
  }

  auto testPow = [](const double& number, const int& exponent) -> bool {
    const double test = ConstexprMagic::Math::pow(number, exponent);
    const double reference = std::pow(number, exponent);

    bool passes = ConstexprMagic::floating::isCloseRelative(
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
    TemplateMagic::all_of(
      TemplateMagic::zipMap(
        TemplateMagic::random.getN<double>(-1e5, 1e5, numTests),
        TemplateMagic::random.getN<int>(-40, 40, numTests),
        testPow
      )
    )
  );

  // ln
  const auto randomZ = TemplateMagic::random.getN<double>(1e-10, 1e10, numTests);
  bool all_ln_pass = TemplateMagic::all_of(
    TemplateMagic::map(
      randomZ,
      [](const auto& z) -> bool {
        bool pass = ConstexprMagic::floating::isCloseRelative(
          ConstexprMagic::Math::ln(z),
          std::log(z),
          accuracy
        );

        if(!pass) {
          std::cout << "ln deviates for z = " << std::setw(12) << z 
            << ", ln(z) = " << std::setw(12) << ConstexprMagic::Math::ln(z) 
            << ", std::log(z) = " << std::setw(12) << std::log(z)
            << ", |Δ| = " << std::setw(12) << std::fabs(
              ConstexprMagic::Math::ln(z) - std::log(z)
            ) << std::endl;
        }

        return pass;
      }
    )
  );

  BOOST_CHECK(all_ln_pass);

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::map(
        TemplateMagic::random.getN<double>(-100, 100, numTests),
        [](const double& x) -> bool {
          return(ConstexprMagic::Math::floor(x) <= x);
        }
      )
    )
  );

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::map(
        TemplateMagic::random.getN<double>(-100, 100, numTests),
        [](const double& x) -> bool {
          return(ConstexprMagic::Math::ceil(x) >= x);
        }
      )
    )
  );

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::map(
        TemplateMagic::random.getN<double>(-100, 100, numTests),
        [](const double& x) -> bool {
          const double rounded = ConstexprMagic::Math::round(x);
          return(
            rounded == ConstexprMagic::Math::floor(x)
            || rounded == ConstexprMagic::Math::ceil(x)
          );
        }
      )
    )
  );

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::map(
        TemplateMagic::random.getN<double>(-M_PI / 2, M_PI / 2, numTests),
        [](const double& x) -> bool {
          return ConstexprMagic::floating::isCloseRelative(
            ConstexprMagic::Math::atan(x),
            std::atan(x),
            accuracy
          );
        }
      )
    )
  );
}

BOOST_AUTO_TEST_CASE(arrayPermutation) {
  std::array<unsigned, 4> base {0, 1, 2, 3};
  std::array<unsigned, 4> STLComparison {0, 1, 2, 3};

  bool customHasNext = true;
  bool STLHasNext = true;

  do {
    customHasNext = ConstexprMagic::inPlaceNextPermutation(base);
    STLHasNext = std::next_permutation(STLComparison.begin(), STLComparison.end());

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In forward permutation, base is {" << TemplateMagic::condenseIterable(base)
      << "} and and STL is {" << TemplateMagic::condenseIterable(STLComparison)
      << "}"
    );
  } while(customHasNext && STLHasNext);

  BOOST_CHECK_MESSAGE(
   !customHasNext && !STLHasNext,
   "The two permutation algorithms don't terminate at the same time"
  );

  base = {3, 2, 1, 0};
  STLComparison = {3, 2, 1, 0};
  customHasNext = true;
  STLHasNext = true;

  do {
    customHasNext = ConstexprMagic::inPlacePreviousPermutation(base);
    STLHasNext = std::prev_permutation(STLComparison.begin(), STLComparison.end());

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In backward permutation, base is {" << TemplateMagic::condenseIterable(base)
      << "} and and STL is {" << TemplateMagic::condenseIterable(STLComparison)
      << "}"
    );
  } while(customHasNext && STLHasNext);

  BOOST_CHECK_MESSAGE(
   !customHasNext && !STLHasNext,
   "The two permutation algorithms don't terminate at the same time"
  );

  // Variants of limited sections of the array

  base = {0, 1, 2, 3};
  STLComparison = {0, 1, 2, 3};
  customHasNext = true;
  STLHasNext = true;

  do {
    customHasNext = ConstexprMagic::inPlaceNextPermutation(base, 1, 3);
    STLHasNext = std::next_permutation(STLComparison.begin() + 1, STLComparison.end() - 1);

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In limited forward permutation, base is {" << TemplateMagic::condenseIterable(base)
      << "} and and STL is {" << TemplateMagic::condenseIterable(STLComparison)
      << "}"
    );
  } while(customHasNext && STLHasNext);

  BOOST_CHECK_MESSAGE(
   !customHasNext && !STLHasNext,
   "The two permutation algorithms don't terminate at the same time in limited "
   "forward permutation"
  );

  base = {3, 2, 1, 0};
  STLComparison = {3, 2, 1, 0};
  customHasNext = true;
  STLHasNext = true;

  do {
    customHasNext = ConstexprMagic::inPlacePreviousPermutation(base, 1, 3);
    STLHasNext = std::prev_permutation(STLComparison.begin() + 1, STLComparison.end() - 1);

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In limited backward permutation, base is {" << TemplateMagic::condenseIterable(base)
      << "} and and STL is {" << TemplateMagic::condenseIterable(STLComparison)
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
  ConstexprMagic::DynamicArray<unsigned, 10> nonConstArr {4, 3, 6};
  nonConstArr.push_back(9);

  return nonConstArr.size() == 4;
}

BOOST_AUTO_TEST_CASE(dynamicArray) {
  constexpr ConstexprMagic::DynamicArray<unsigned, 10> arr {4, 3, 5};

  static_assert(
    arr.size() == 3,
    "Array size isn't initialized correctly from parameter pack ctor"
  );
  static_assert(
    compileTimeDynTest(),
    "non-const dynamic array functionality works as expected"
  );

  constexpr ConstexprMagic::Array<unsigned, 10> values {1, 2, 2, 3, 3, 3, 4, 4, 4, 4};

  constexpr auto grouped = groupByEquality(
    values,
    std::equal_to<unsigned>()
  );

  static_assert(
    grouped.size() == 4
    && grouped.at(0).size() == 1
    && grouped.at(1).size() == 2
    && grouped.at(2).size() == 3
    && grouped.at(3).size() == 4,
    "Grouping does not work as expected"
  );
}

template<typename T, size_t size, class Comparator>
constexpr bool isSorted(const ConstexprMagic::DynamicSet<T, size, Comparator>& set) {
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

BOOST_AUTO_TEST_CASE(dynamicSet) {
  ConstexprMagic::DynamicSet<unsigned, 10> set;

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
        << TemplateMagic::condenseIterable(set)
        << "}."
    );
  }
  BOOST_CHECK(isSorted(set));

  ConstexprMagic::DynamicSet<unsigned, 10> setInitList {
    ConstexprMagic::DynamicArray<unsigned, 10> {
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
      << TemplateMagic::condenseIterable(setInitList)
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
  ConstexprMagic::Array<
    ConstexprMagic::Array<unsigned, 4>,
    5
  > sampleArrays {
    ConstexprMagic::Array<unsigned, 4> {1u, 2u, 3u, 4u},
    ConstexprMagic::Array<unsigned, 4> {1u, 2u, 4u, 3u},
    ConstexprMagic::Array<unsigned, 4> {1u, 4u, 3u, 2u},
    ConstexprMagic::Array<unsigned, 4> {1u, 4u, 2u, 3u},
    ConstexprMagic::Array<unsigned, 4> {2u, 1u, 3u, 4u}
  };

  ConstexprMagic::DynamicSet<
    ConstexprMagic::Array<unsigned, 4>,
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
bool validate(const ConstexprMagic::DynamicSet<T, size>& set) {
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
    set.size() == std::distance(
      set.begin(),
      set.end()
    )
  );
}

BOOST_AUTO_TEST_CASE(arrayOperators) {
  ConstexprMagic::Array<unsigned, 4> a {4, 2, 3, 1};
  ConstexprMagic::Array<unsigned, 4> b {4, 3, 2, 1};
  
  BOOST_CHECK(ConstexprMagic::testLogicalOperators(a, b));
  BOOST_CHECK(ConstexprMagic::testLogicalOperators(a, a));

  ConstexprMagic::dynamic::explainLogicalOperatorFailures(a, b);
}

BOOST_AUTO_TEST_CASE(dynamicSetFuzzing) {
  for(unsigned N = 0; N < 100; ++N) {
    ConstexprMagic::DynamicSet<unsigned, 100> subject;

    std::vector<unsigned> numbers;
    numbers.resize(50);
    std::iota(
      numbers.begin(),
      numbers.end(),
      0
    );
    std::shuffle(
      numbers.begin(),
      numbers.end(),
      TemplateMagic::random.randomEngine
    );

    for(const auto& number : numbers) {
      subject.insert(number);

      bool isValid = validate(subject);
      if(!isValid) {
        std::cout << "After inserting " << number 
          << ", set is left in invalid state. Set: {"
          << TemplateMagic::condenseIterable(subject) << "}\ninsert sequence {"
          << TemplateMagic::condenseIterable(numbers) << "}."
          << std::endl;
      }

      BOOST_REQUIRE(isValid);

      BOOST_CHECK(subject.contains(number));
    }
  }
}

BOOST_AUTO_TEST_CASE(specificTest) {
  ConstexprMagic::DynamicSet<
    ConstexprMagic::Array<unsigned, 6>,
    5
  > set {
    ConstexprMagic::Array<unsigned, 6> {0, 1, 2, 3, 4, 5}
  };

  BOOST_CHECK(
    !set.contains(
      ConstexprMagic::Array<unsigned, 6> {1, 2, 3, 4, 0, 5}
    )
  );
}

namespace TupleTypeTests {

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
  ConstexprMagic::TupleType::unpackToFunction<Fruit, sumNumbersValue>() == 57,
  "Unpacking fruit tuple to valueValue does not yield expected result"
);

static_assert(
  ConstexprMagic::TupleType::unpackToFunction<Fruit, sumNumbersFunctor>() == 57,
  "Unpacking fruit tuple to valueFunctor does not yield expected result"
);

static_assert(
  ConstexprMagic::arraysEqual(
    ConstexprMagic::TupleType::map<Fruit, getNumberValue>(),
    std::array<unsigned, 3> {{Apple::number, Banana::number, Cherry::number }}
  ),
  "Mapping with getNumberValue does not yield expected result!"
);

static_assert(
  ConstexprMagic::arraysEqual(
    ConstexprMagic::TupleType::map<Fruit, getNumberFunctor>(),
    std::array<unsigned, 3> {{Apple::number, Banana::number, Cherry::number }}
  ),
  "Mapping with getNumberFunctor does not yield expected result!"
);

static_assert(
  ConstexprMagic::arraysEqual(
    ConstexprMagic::TupleType::mapAllPairs<Fruit, pairSumValue>(),
    std::array<unsigned, 3> {{
      Apple::number + Banana::number,
      Apple::number + Cherry::number,
      Banana::number + Cherry::number
    }}
  ),
  "Mapping pairs with pairSumValue does not yield expected result!"
);

static_assert(
  ConstexprMagic::arraysEqual(
    ConstexprMagic::TupleType::mapAllPairs<Fruit, pairSumFunctor>(),
    std::array<unsigned, 3> {{
      Apple::number + Banana::number,
      Apple::number + Cherry::number,
      Banana::number + Cherry::number
    }}
  ),
  "Mapping pairs with pairSumFunctor does not yield expected result!"
);

using countTestType = std::tuple<unsigned, float, double, unsigned, size_t>;

static_assert(
  ConstexprMagic::TupleType::countType<countTestType, unsigned>() == 2,
  "Counting unsigned in countTestType does not return two!"
);

} // namespace TupleTypeTests

namespace FloatingPointComparisonTests {

template<typename T>
constexpr bool testAbsoluteComparison(const T& a, const T& b, const T& tolerance) {
  ConstexprMagic::floating::ExpandedAbsoluteEqualityComparator<T> comparator {
    tolerance
  };

  return (
    ConstexprMagic::Math::XOR(
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
    ) && ConstexprMagic::Math::XOR(
      comparator.isEqual(a, b),
      comparator.isUnequal(a, b)
    )
  );
}

template<typename T>
constexpr bool testRelativeComparison(const T& a, const T& b, const T& tolerance) {
  ConstexprMagic::floating::ExpandedRelativeEqualityComparator<T> comparator {
    tolerance
  };

  return (
    ConstexprMagic::Math::XOR(
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
    ) && ConstexprMagic::Math::XOR(
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

constexpr ConstexprMagic::Array<unsigned, 4> f {4, 2, 9, 3};
constexpr ConstexprMagic::Array<unsigned, 4> g {11, 22, 33, 44};
constexpr ConstexprMagic::Array<unsigned, 4> h {234, 292, 912, 304};
constexpr ConstexprMagic::Array<unsigned, 8> fg {
  4, 2, 9, 3,
  11, 22, 33, 44
};
constexpr ConstexprMagic::Array<unsigned, 12> fgh {
  4, 2, 9, 3,
  11, 22, 33, 44,
  234, 292, 912, 304
};

static_assert(
  ConstexprMagic::arraysEqual(
    ConstexprMagic::arrayConcatenate(f, g),
    fg
  ),
  "Pairwise concatenation does not preserve sequence!"
);

static_assert(
  ConstexprMagic::arraysEqual(
    ConstexprMagic::arrayConcatenate(f, g, h),
    fgh
  ),
  "Variadic concatenation does not work as expected"
);

} // namespace ConcatenationTests
