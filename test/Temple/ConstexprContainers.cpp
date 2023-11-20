/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"

#include "Molassembler/Temple/constexpr/Array.h"
#include "Molassembler/Temple/constexpr/Bitmask.h"
#include "Molassembler/Temple/constexpr/Bitset.h"
#include "Molassembler/Temple/constexpr/Containers.h"
#include "Molassembler/Temple/constexpr/DynamicArray.h"
#include "Molassembler/Temple/constexpr/DynamicMap.h"
#include "Molassembler/Temple/constexpr/DynamicSet.h"
#include "Molassembler/Temple/constexpr/Jsf.h"
#include "Molassembler/Temple/constexpr/LogicalOperatorTests.h"
#include "Molassembler/Temple/constexpr/Math.h"
#include "Molassembler/Temple/constexpr/TupleType.h"
#include "Molassembler/Temple/constexpr/TupleTypePairs.h"
#include "Molassembler/Temple/Permutations.h"

#include <iostream>

using namespace Scine::Molassembler;
extern Temple::Generator<> generator;

using namespace std::string_literals;

namespace ArrayTests {

constexpr auto testArr = Temple::Array<unsigned, 3> {4, 3, 5};

template<typename T, size_t size>
constexpr Temple::Array<T, size> modifyArray(
  const Temple::Array<T, size>& array
) {
  auto arrayCopy = array;
  inPlaceSwap(arrayCopy, 0, 1);
  return arrayCopy;
}

constexpr auto modf = modifyArray(testArr);

static_assert(
  modf == Temple::Array<unsigned, 3> {3, 4, 5},
  "Swap doesn't work as expected"
);

static_assert(
  std::is_same<
    decltype(Temple::makeArray(4, 3, 9)),
    Temple::Array<int, 3>
  >::value,
  "makeArray does not work as expected"
);

} // namespace ArrayTests

BOOST_AUTO_TEST_CASE(ArrayPermutation, *boost::unit_test::label("Temple")) {
  std::array<unsigned, 4> base {{0, 1, 2, 3}};
  std::array<unsigned, 4> STLComparison {{0, 1, 2, 3}};

  bool customHasNext;
  bool STLHasNext;

  do {
    customHasNext = Temple::inPlaceNextPermutation(base);
    STLHasNext = std::next_permutation(STLComparison.begin(), STLComparison.end());

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In forward permutation, base is {" << Temple::condense(base)
      << "} and and STL is {" << Temple::condense(STLComparison)
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
    customHasNext = Temple::inPlacePreviousPermutation(base);
    STLHasNext = std::prev_permutation(STLComparison.begin(), STLComparison.end());

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In backward permutation, base is {" << Temple::condense(base)
      << "} and and STL is {" << Temple::condense(STLComparison)
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
    customHasNext = Temple::inPlaceNextPermutation(base, 1, 3);
    STLHasNext = std::next_permutation(STLComparison.begin() + 1, STLComparison.end() - 1);

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In limited forward permutation, base is {" << Temple::condense(base)
      << "} and and STL is {" << Temple::condense(STLComparison)
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
    customHasNext = Temple::inPlacePreviousPermutation(base, 1, 3);
    STLHasNext = std::prev_permutation(STLComparison.begin() + 1, STLComparison.end() - 1);

    BOOST_CHECK_MESSAGE(
      base == STLComparison,
      "In limited backward permutation, base is {" << Temple::condense(base)
      << "} and and STL is {" << Temple::condense(STLComparison)
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
  Temple::DynamicArray<unsigned, 10> nonConstArr {4, 3, 6};
  nonConstArr.push_back(9);
  return nonConstArr.size() == 4;
}

constexpr bool dynArrSpliceTest() {
  Temple::DynamicArray<unsigned, 10> nonConstArr {4, 3, 6, 5, 1, 9};
  auto spliced = nonConstArr.splice(2);

  return (
    spliced == Temple::DynamicArray<unsigned, 10> {6, 5, 1, 9}
    && nonConstArr == Temple::DynamicArray<unsigned, 10> {4, 3}
  );
}

BOOST_AUTO_TEST_CASE(DynamicArrayTests, *boost::unit_test::label("Temple")) {
  constexpr Temple::DynamicArray<unsigned, 10> arr {4, 3, 5};

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

  constexpr Temple::Array<unsigned, 10> values {1, 2, 2, 3, 3, 3, 4, 4, 4, 4};

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

  constexpr Temple::DynamicArray<unsigned, 14> fromFixed {values};

  static_assert(fromFixed.size() == 10, "Construction from fixed doesn't work");
}

template<typename T, size_t size, class Comparator>
constexpr bool isSorted(const Temple::DynamicSet<T, size, Comparator>& set) {
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

BOOST_AUTO_TEST_CASE(DynamicSetTests, *boost::unit_test::label("Temple")) {
  Temple::DynamicSet<unsigned, 10> set;

  BOOST_CHECK(set.size() == 0);
  BOOST_CHECK(
    std::distance(
      set.begin(),
      set.end()
    ) == 0
  );

  for(const auto& item : {9u, 3u, 5u}) {
    set.insert(item);
  }

  BOOST_CHECK(set.size() == 3);
  BOOST_CHECK(
    std::distance(
      set.begin(),
      set.end()
    ) == 3
  );
  BOOST_CHECK(set.contains(3) && set.contains(5) && set.contains(9));
  for(const auto& item : {2u, 4u, 8u, 10u}) {
    BOOST_CHECK_MESSAGE(
      !set.contains(item),
      "Set says it contains " << item << " when it shouldn't (set is {"
        << Temple::condense(set)
        << "}."
    );
  }
  BOOST_CHECK(isSorted(set));

  Temple::DynamicSet<unsigned, 10> setInitList {
    Temple::DynamicArray<unsigned, 10> {4, 9, 13}
  };

  BOOST_CHECK(setInitList.size() == 3);
  BOOST_CHECK(
    std::distance(
      setInitList.begin(),
      setInitList.end()
    ) == 3
  );

  setInitList.insert(0);

  BOOST_CHECK(setInitList.size() == 4);

  BOOST_CHECK_MESSAGE(
    setInitList.contains(4)
    && setInitList.contains(9)
    && setInitList.contains(13)
    && setInitList.contains(0)
    && !setInitList.contains(1)
    && !setInitList.contains(25),
    "setInitList {"
      << Temple::condense(setInitList)
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
  Temple::Array<
    Temple::Array<unsigned, 4>,
    5
  > sampleArrays {
    Temple::Array<unsigned, 4> {1, 2, 3, 4},
    Temple::Array<unsigned, 4> {1, 2, 4, 3},
    Temple::Array<unsigned, 4> {1, 4, 3, 2},
    Temple::Array<unsigned, 4> {1, 4, 2, 3},
    Temple::Array<unsigned, 4> {2, 1, 3, 4}
  };

  Temple::DynamicSet<
    Temple::Array<unsigned, 4>,
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
bool validate(const Temple::DynamicSet<T, size>& set) {
  /* NOTE: std::is_sorted isn't constexpr yet, and the range adaptors
   * aren't constexpr, so we need to do the sequential pairs thing ourselves
   */

  // Is the set ordered?
  auto leftIter = set.begin();
  auto rightIter = leftIter; ++rightIter;
  auto bound = set.end(); --bound;

  if(leftIter == set.end() || rightIter == set.end()) {
    return true;
  }

  while(rightIter != bound) {
    BOOST_REQUIRE_LT(*leftIter, *rightIter);

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

BOOST_AUTO_TEST_CASE(ArrayOperators, *boost::unit_test::label("Temple")) {
  Temple::Array<unsigned, 4> a {4, 2, 3, 1};
  Temple::Array<unsigned, 4> b {4, 3, 2, 1};

  BOOST_CHECK(Temple::testLogicalOperators(a, b));
  BOOST_CHECK(Temple::testLogicalOperators(a, a));

  Temple::dynamic::explainLogicalOperatorFailures(a, b);
}

BOOST_AUTO_TEST_CASE(DynamicSetFuzzing, *boost::unit_test::label("Temple")) {
  for(unsigned N = 0; N < 100; ++N) {
    Temple::DynamicSet<unsigned, 100> subject;

    std::vector<unsigned> numbers;
    numbers.resize(50);
    std::iota(
      numbers.begin(),
      numbers.end(),
      0
    );

    Temple::Random::shuffle(numbers, generator.engine);

    for(const auto& number : numbers) {
      subject.insert(number);

      const bool isValid = validate(subject);
      if(!isValid) {
        std::cout << "After inserting " << number
          << ", set is left in invalid state. Set: {"
          << Temple::condense(subject) << "}\ninsert sequence {"
          << Temple::condense(numbers) << "}."
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
  Temple::Tuples::unpackToFunction<Fruit, sumNumbersValue>() == 57,
  "Unpacking fruit tuple to valueValue does not yield expected result"
);

static_assert(
  Temple::Tuples::unpackToFunction<Fruit, sumNumbersFunctor>() == 57,
  "Unpacking fruit tuple to valueFunctor does not yield expected result"
);

static_assert(
  Temple::arraysEqual(
    Temple::Tuples::map<Fruit, getNumberValue>(),
    Temple::Array<unsigned, 3> {{Apple::number, Banana::number, Cherry::number}}
  ),
  "Mapping with getNumberValue does not yield expected result!"
);

static_assert(
  Temple::arraysEqual(
    Temple::Tuples::map<Fruit, getNumberFunctor>(),
    Temple::Array<unsigned, 3> {{Apple::number, Banana::number, Cherry::number}}
  ),
  "Mapping with getNumberFunctor does not yield expected result!"
);

static_assert(
  Temple::arraysEqual(
    Temple::Tuples::mapAllPairs<Fruit, pairSumValue>(),
    Temple::Array<unsigned, 3> {{
      Apple::number + Banana::number,
      Apple::number + Cherry::number,
      Banana::number + Cherry::number
    }}
  ),
  "Mapping pairs with pairSumValue does not yield expected result!"
);

static_assert(
  Temple::arraysEqual(
    Temple::Tuples::mapAllPairs<Fruit, pairSumFunctor>(),
    Temple::Array<unsigned, 3> {{
      Apple::number + Banana::number,
      Apple::number + Cherry::number,
      Banana::number + Cherry::number
    }}
  ),
  "Mapping pairs with pairSumFunctor does not yield expected result!"
);

using countTestType = std::tuple<unsigned, float, double, unsigned, size_t>;

static_assert(
  Temple::Tuples::countType<countTestType, unsigned>() == 2,
  "Counting unsigned in countTestType does not return two!"
);

} // namespace tuple_type_tests

namespace ConcatenationTests {

constexpr Temple::Array<unsigned, 4> f {4, 2, 9, 3};
constexpr Temple::Array<unsigned, 4> g {11, 22, 33, 44};
constexpr Temple::Array<unsigned, 4> h {234, 292, 912, 304};
constexpr Temple::Array<unsigned, 8> fg {
  4, 2, 9, 3,
  11, 22, 33, 44
};
constexpr Temple::Array<unsigned, 12> fgh {
  4, 2, 9, 3,
  11, 22, 33, 44,
  234, 292, 912, 304
};

static_assert(
  Temple::arraysEqual(
    Temple::arrayConcatenate(f, g),
    fg
  ),
  "Pairwise concatenation does not preserve sequence!"
);

static_assert(
  Temple::arraysEqual(
    Temple::arrayConcatenate(f, g, h),
    fgh
  ),
  "Variadic concatenation does not work as expected"
);

} // namespace ConcatenationTests

namespace DynamicMapTests {

constexpr Temple::DynamicMap<unsigned, int, 20> generateMap() {
  Temple::DynamicMap<unsigned, int, 20> myMap;

  myMap.insert(4, -2);
  myMap.insert(1, 4);
  myMap.insert(3, 9);

  return myMap;
}

constexpr auto a = generateMap();

static_assert(a.at(4) == -2, "Map does not find element with key 4");
static_assert(a.at(1) == 4, "Map does not find element with key 1");
static_assert(a.at(3) == 9, "Map does not find element with key 3");

} // namespace DynamicMapTests

namespace UpperTriangularMatrixTests {

// Can default-construct
constexpr auto defaultMatr = Temple::UpperTriangularMatrix<bool, 15> {};
static_assert(decltype(defaultMatr)::N == 6, "Size isn't right");

constexpr auto matr = Temple::makeUpperTriangularMatrix(
  std::array<unsigned, 6> {{1, 2, 3, 4, 5, 6}}
);

static_assert(decltype(matr)::N == 4, "Size isn't right");

/*constexpr auto failing = Temple::makeUpperTriangularMatrix(
  std::array<unsigned, 5> {{1, 2, 3, 4, 5}}
);*/

constexpr auto fromArray = Temple::makeUpperTriangularMatrix(
  Temple::Array<unsigned, 6> {{1, 2, 3, 4, 5, 6}}
);
static_assert(decltype(fromArray)::N == 4, "Size isn't right");

} // namespace UpperTriangularMatrixTests

enum class ScopedEnum : unsigned {A, B, C};
enum UnscopedEnum : unsigned {D, E};

BOOST_AUTO_TEST_CASE(BitmaskAll, *boost::unit_test::label("Temple")) {
  using namespace Temple;

  constexpr auto a = make_bitmask(ScopedEnum::A) | ScopedEnum::C;
  static_assert(a & ScopedEnum::A, "A must be contained in the bitmask a");
  BOOST_CHECK(a[ScopedEnum::A] && a[ScopedEnum::C] && !a[ScopedEnum::B]);

  constexpr auto b = make_bitmask(UnscopedEnum::E);
  BOOST_CHECK(b[UnscopedEnum::E] && !b[UnscopedEnum::D]);
}

constexpr Temple::Bitset<300> make_a_bitset() {
  Temple::Bitset<300> a;

  a.set(0);
  a.set(299);
  a.set(63);
  a.set(128);

  return a;
}

BOOST_AUTO_TEST_CASE(BitsetTests, *boost::unit_test::label("Temple")) {
  using namespace Temple;

  constexpr auto a = make_a_bitset();

  static_assert(a.test(0), "Zero hasn't been set");

  BOOST_CHECK(a.test(0) && a.test(299) && a.test(63) && a.test(128));
  BOOST_CHECK(!a.test(1) && !a.test(298) && !a.test(64) && !a.test(127));

  auto b = make_a_bitset();
  b.unset(63);
  BOOST_CHECK(!b.test(63));
}

BOOST_AUTO_TEST_CASE(PermutationIndexTests, *boost::unit_test::label("Temple")) {
  using namespace Temple;

  auto a = iota<Array, unsigned, 6>();

  // expect monotonically increasing permutation index on inplaceNextPermutation
  auto index = permutationIndex(a);

  BOOST_CHECK_EQUAL(index, 0);

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
