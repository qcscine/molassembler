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
      return ConstexprMagic::Math::isCloseRelative(
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
      return ConstexprMagic::Math::isCloseRelative(
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

    bool passes = ConstexprMagic::Math::isCloseRelative(
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
        bool pass = ConstexprMagic::Math::isCloseRelative(
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

  for(const auto& item : {9u, 3u, 5u}) {
    set.insert(item);
  }

  BOOST_CHECK(set.size() == 3);
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
}
