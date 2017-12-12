#define BOOST_TEST_MODULE ConnectivityManagerTests
#define BOOST_TEST_DYN_LINK

#include <boost/test/results_collector.hpp>
#include <boost/test/unit_test.hpp>

#include "BTree.h"
#include "Containers.h"
#include "Enumerate.h"
#include "MemberFetcher.h"
#include "Numeric.h"
#include "Random.h"
#include "VectorView.h"
#include "Optionals.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

template<typename T>
std::ostream& operator << (std::ostream& os, const boost::optional<T>& valueOptional) {
  if(valueOptional) {
    os << "Some " << valueOptional.value();
  } else {
    os << "None";
  }

  return os;
}

inline bool lastTestPassed() {
  using namespace boost::unit_test;

  test_case::id_t id = framework::current_test_case().p_id;
  test_results rez = results_collector.results(id);
  return rez.passed();
}

template<typename NullaryCallable>
double timeNullaryCallable(
  NullaryCallable&& function
) {
  using namespace std::chrono;

  time_point<system_clock> start, end;
  std::array<unsigned, 1000> timings;
  for(unsigned N = 0; N < 1000; ++N) {
    start = system_clock::now();

    function();

    end = system_clock::now();
    duration<double> elapsed = end - start;

    timings.at(N) = elapsed.count() * 1e9;
  }

  return TemplateMagic::average(timings);
}

double divByThree (unsigned a) {
  return static_cast<double>(a) / 3.0;
}

BOOST_AUTO_TEST_CASE( sumTest ) {
  std::vector<unsigned> instance {0, 1, 2, 3};
  auto f = TemplateMagic::sum(instance);

  BOOST_CHECK(f == 6);

  auto mapped = TemplateMagic::map(
    instance,
    divByThree
  );

  BOOST_CHECK(mapped == std::vector<double>({0, 1.0/3.0, 2.0/3.0, 1}));

  auto pairwiseSum = TemplateMagic::mapSequentialPairs(
    instance,
    std::plus<unsigned>()
  );

  BOOST_CHECK(pairwiseSum == std::vector<unsigned>({1,3,5}));

  auto pairwiseSmaller = TemplateMagic::accumulate(
    TemplateMagic::mapSequentialPairs(
      instance,
      std::less<unsigned>()
    ),
    true,
    std::logical_and<bool>()
  );

  BOOST_CHECK(pairwiseSmaller);

  std::vector<
    std::vector<unsigned>
  > vectorOfVectors {
    {0, 1, 4},
    {4, 5}
  };

  auto mapToSizes = TemplateMagic::map(
    vectorOfVectors,
    [](const std::vector<unsigned>& vectorUnsigned) -> unsigned {
      return vectorUnsigned.size();
    }
  );

  std::vector<unsigned> unsignedVector {1, 2, 3};

  BOOST_CHECK(
    TemplateMagic::sum(
      TemplateMagic::mapAllPairs(
        unsignedVector,
        [](const unsigned& a, const unsigned& b) -> unsigned {
          return a + b;
        }
      )
    ) == 12
  );

  std::vector<double> doubleVector {1.2, 1.5, 1.9};

  BOOST_CHECK(
    TemplateMagic::sum(
      TemplateMagic::mapAllPairs(
        doubleVector,
        [](const double& a, const double& b) -> double {
          return a + b;
        }
      )
    ) == 9.2
  );
}

BOOST_AUTO_TEST_CASE( reduceTests) {
  std::vector<unsigned> values {1, 2, 3, 4, 5};
  BOOST_CHECK(
    TemplateMagic::reduce(
      values,
      0u,
      std::plus<unsigned>()
    ) == 15u
  );
  BOOST_CHECK(
    TemplateMagic::reduce(
      values,
      1u,
      std::multiplies<unsigned>()
    ) == 120u
  );
}

BOOST_AUTO_TEST_CASE( minMaxTests ) {
  const std::vector<unsigned> values {1, 4, 6, 8};
  BOOST_CHECK(TemplateMagic::max(values) == 8u);
  BOOST_CHECK(TemplateMagic::min(values) == 1u);
}

BOOST_AUTO_TEST_CASE(kahanSummation) {
  for(unsigned nTest = 0; nTest < 100; nTest++) {
    const unsigned N = 100;
    const unsigned magnitudeSpread = 20;
    const auto randomNumbers = TemplateMagic::random.getN<double>(
      std::pow(10, - static_cast<double>(magnitudeSpread) / 2),
      std::pow(10, static_cast<double>(magnitudeSpread) / 2),
      N
    );

    const double reduceSum = TemplateMagic::reduce(
      randomNumbers,
      0.0,
      std::plus<double>()
    );

    const double kahanSum = TemplateMagic::kahanSum(randomNumbers);

    /* Reference sum with long doubles, I know an alternative implementation of
     * standard reduce summation with long intermediates would have done the
     * trick too but I'm lazy and this is just a test
     */
    const auto addedPrecision = TemplateMagic::cast<long double>(randomNumbers);
    const double longSum = TemplateMagic::reduce(
      addedPrecision,
      0.0l,
      std::plus<long double>()
    );

    // Kahan summation should be equally or more accurate than the standard reduction
    BOOST_CHECK_MESSAGE(
      std::fabs(longSum - reduceSum) >= std::fabs(longSum - kahanSum),
      "Kahan summation is less accurate than standard reduce sum! "
      << "long: " << longSum << ", kahan: " << kahanSum 
      << ", reduce: " << reduceSum << ". Absolute deviation from long sum: "
      << "kahan:" << std::fabs(longSum - kahanSum) << ", "
      << "reduce: " << std::fabs(longSum - reduceSum)
    );
  }
}

BOOST_AUTO_TEST_CASE(numericAverageStdDev) {
  const std::vector<double> values {29, 30, 31, 32, 33};

  BOOST_CHECK(TemplateMagic::average(values) == 31); 
  BOOST_CHECK(
    std::fabs(TemplateMagic::stddev(values) - std::sqrt(2)) 
    < 1e-10
  );
}

BOOST_AUTO_TEST_CASE(mapToSameContainerTests) {
  std::set<int> f {5, -1, 9};

  auto fMapped = TemplateMagic::map(
    f,
    [](const int& x) -> double {
      return x + 1.3;
    }
  );

  static_assert(
    std::is_same<decltype(fMapped), std::set<double>>::value, 
    "Map to same container does not work as expected"
  );

  std::vector<float> x {0, 3.4, 9};
  auto xMapped = TemplateMagic::map(
    x,
    [](const float& x) -> unsigned long {
      return static_cast<unsigned long>(x);
    }
  );

  static_assert(
    std::is_same<decltype(xMapped), std::vector<unsigned long>>::value, 
    "Map to same container does not work as expected"
  );
}

BOOST_AUTO_TEST_CASE( boolArrayAllOf) {
  std::array<bool, 3> testArray {true, false, true};

  BOOST_CHECK(!TemplateMagic::all_of(testArray));
}

/*BOOST_AUTO_TEST_CASE( enumerateTests) {
  std::vector<unsigned> testVec {5, 2, 3, 4};

  std::cout << "Before enumerate:" << std::endl;
  for(const auto& enumStruct : enumerate(testVec)) {
    std::cout << "{ index: " << enumStruct.index << ", value: " 
      << enumStruct.value << "}" << std::endl;
  }
}*/

BOOST_AUTO_TEST_CASE(memberFetcherTests) {
  std::vector<
    std::vector<unsigned>
  > testClass {
    {4, 3, 9, 10},
    {1, 2, 3, 4, 5},
    {11},
    {44, 93}
  };

  auto memberFetcherClass = TemplateMagic::getMember(
    testClass,
    [](const auto& container) -> unsigned {
      return container.size();
    }
  );

  BOOST_CHECK_MESSAGE(
    TemplateMagic::sum(memberFetcherClass) == 12u,
    "Member fetcher does not work as expected, begin()-end() contains {"
      << TemplateMagic::condenseIterable(memberFetcherClass)
      << "} instead of expected {4, 5, 1, 2}"
  );

  BOOST_CHECK_MESSAGE(
    TemplateMagic::min(memberFetcherClass) == 1
    && TemplateMagic::max(memberFetcherClass) == 5,
    "Usage of min/max does not return expected result"
  );

  // Is this approach more efficient than copying member information?
  double copyingTime = timeNullaryCallable(
    [&]() {
      auto mappedData = TemplateMagic::map(
        testClass,
        [](const auto& subVector) -> unsigned {
          return subVector.size();
        }
      );

      for(const auto& size : mappedData) {
        static_cast<void>(size);
      }
    }
  );

  double memberFetcherTime = timeNullaryCallable(
    [&]() {
      auto fetcher = TemplateMagic::getMember(
        testClass,
        [](const auto& subVector) -> unsigned {
          return subVector.size();
        }
      );

      for(const auto& size : fetcher) {
        static_cast<void>(size);
      }
    }
  );

  BOOST_CHECK_MESSAGE(
    memberFetcherTime < copyingTime,
    "Copying is faster than fetching members: copy = " << copyingTime << " vs. " 
    << memberFetcherTime << " member fetch " 
  );

  // Interoperability with VectorView

  auto filteredView = TemplateMagic::filter(
    testClass,
    [](const auto& subVector) -> bool {
      return subVector.size() < 3;
    }
  );

  auto minSize = TemplateMagic::min(
    TemplateMagic::getMember(
      filteredView,
      [](const auto& subVector) -> unsigned {
        return subVector.size();
      }
    )
  );

  BOOST_CHECK_MESSAGE(
    minSize == 4,
    "Min size of filteredView returns " << minSize << ", not 4"
  );

  auto maxSize = TemplateMagic::max(
    TemplateMagic::getMember(
      filteredView,
      [](const auto& subVector) -> unsigned {
        return subVector.size();
      }
    )
  );

  BOOST_CHECK_MESSAGE(
    maxSize == 5,
    "Max size of filteredView returns " << maxSize << ", not 5"
  );
}

BOOST_AUTO_TEST_CASE(concatenateTests) {
  std::vector<unsigned> f {{4, 9, 1}};
  std::set<unsigned> x {{3, 6, 10}};

  std::vector<unsigned> comp {{4, 9, 1, 3, 6, 10}};

  static_assert(
    std::is_same<
      decltype(TemplateMagic::concatenate(f, x)),
      std::vector<unsigned>
    >::value,
    "Concatenate must return a vector of the underlying shared type"
  );

  BOOST_CHECK(
    TemplateMagic::concatenate(f, x)
    == comp
  );
}

namespace BTreeStaticTests {

TemplateMagic::BTree<unsigned, 3> generateTree() {
  TemplateMagic::BTree<unsigned, 3> tree;

  tree.insert(9);
  tree.insert(3);
  tree.insert(5);
  tree.insert(20);

  return tree;
}

auto testTree = generateTree();

static_assert(
  /* BTree of minimum order 3 has max 5 keys per node and max 6 children per node
   *
   * height  nodes       keys
   * 0       1           5
   * 1       1 + 6       5 + 6*5
   * 2       1 + 6 + 36  5 + 6*5 + 36*5
   *
   * #nodes(h) = sum_{i = 0}^{h} (2t)^i
   *
   *     (2t)^{h + 1} - 1
   *  N = ----------------
   *         2t - 1
   *
   * -> N * (2t - 1) + 1 = (2t)^{h + 1}
   *
   * -> log_2t [N * (2t - 1) + 1] = h + 1
   *
   * -> h = log_2t [N * (2t - 1) + 1] - 1
   *
   */
  TemplateMagic::BTreeProperties::minHeight(5, 3) == 0
  && TemplateMagic::BTreeProperties::minHeight(35, 3) == 1
  && TemplateMagic::BTreeProperties::minHeight(215, 3) == 2,
  "minHeight function is wrong"
);

static_assert(
  TemplateMagic::BTreeProperties::maxNodesInTree(0, 3) == 1
  && TemplateMagic::BTreeProperties::maxNodesInTree(1, 3) == 7
  && TemplateMagic::BTreeProperties::maxNodesInTree(2, 3) == 43
  && TemplateMagic::BTreeProperties::maxNodesInTree(3, 3) == 259,
  "maxNodesInTree is wrong"
);

} // namespace BTreeStaticTests

unsigned popRandom(std::set<unsigned>& values) {
  auto it = values.begin();

  std::advance(
    it,
    TemplateMagic::random.getSingle<unsigned>(0, values.size() - 1)
  );

  auto value = *it;

  values.erase(it);

  return value;
}

BOOST_AUTO_TEST_CASE(BTreeTests) {
  constexpr unsigned nKeys = 100;

  using namespace std::string_literals;

  std::vector<unsigned> values (nKeys);

  std::iota(
    values.begin(),
    values.end(),
    0
  );

  std::set<unsigned> notInTree {values.begin(), values.end()};
  std::set<unsigned> inTree;

  TemplateMagic::BTree<unsigned, 3> tree;

  std::string lastTreeGraph;

  std::vector<std::string> decisions;

  auto addElement = [&](const std::string& lastTreeGraph) {
    // Add an element
    auto toAdd = popRandom(notInTree);
    decisions.emplace_back("i"s + std::to_string(toAdd));

    BOOST_CHECK_NO_THROW(tree.insert(toAdd));
    BOOST_REQUIRE_MESSAGE(
      lastTestPassed(),
      "Element insertion failed. Operation sequence: "
        << TemplateMagic::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << lastTreeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    inTree.insert(toAdd);
  };

  auto removeElement = [&](const std::string& lastTreeGraph) {
    // Remove an element
    auto toRemove = popRandom(inTree);
    decisions.emplace_back("r"s + std::to_string(toRemove));

    BOOST_CHECK_NO_THROW(tree.remove(toRemove));
    BOOST_REQUIRE_MESSAGE(
      lastTestPassed(),
      "Tree element removal failed. Operation sequence: "
        << TemplateMagic::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << lastTreeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    notInTree.insert(toRemove);
  };

  auto fullValidation = [&](const std::string& lastTreeGraph) {
    // Validate the tree
    BOOST_CHECK_NO_THROW(tree.validate());
    BOOST_REQUIRE_MESSAGE(
      lastTestPassed(),
      "Tree validation failed. Operation sequence: " 
        << TemplateMagic::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << lastTreeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    auto notInsertedNotContained = TemplateMagic::mapToVector(
      notInTree,
      [&](const auto& notInTreeValue) -> bool {
        return !tree.contains(notInTreeValue);
      }
    );

    // Check that all elements are truly contained or not
    BOOST_REQUIRE_MESSAGE(
      TemplateMagic::all_of(notInsertedNotContained),
      "Not all elements recorded as not in the tree are recognized as such!\n" 
        << "Found in the tree, but should not be present: "
        << TemplateMagic::condenseIterable(
          TemplateMagic::moveIf(
            TemplateMagic::zipMap(
              notInsertedNotContained,
              notInTree,
              [](const bool& passed, const unsigned& value) -> std::string {
                if(!passed) {
                  return std::to_string(value);
                }

                return "";
              }
            ),
            [](const std::string& str) -> bool {
              return str != "";
            }
          )
        ) << "\nSequence of operations: " 
        << TemplateMagic::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << lastTreeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    auto insertedContained = TemplateMagic::mapToVector(
      inTree,
      [&](const auto& inTreeValue) -> bool {
        return tree.contains(inTreeValue);
      }
    );

    BOOST_REQUIRE_MESSAGE(
      TemplateMagic::all_of(insertedContained),
      "Not all elements recorded as contained in the tree are recognized as such!\n" 
        << "Not found in the tree: "
        << TemplateMagic::condenseIterable(
          TemplateMagic::moveIf(
            TemplateMagic::zipMap(
              insertedContained,
              inTree,
              [](const bool& passed, const unsigned& value) -> std::string {
                if(!passed) {
                  return std::to_string(value);
                }

                return "";
              }
            ),
            [](const std::string& str) -> bool {
              return str != "";
            }
          )
        ) << "\nSequence of operations: " 
        << TemplateMagic::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << lastTreeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );
  };

  for(unsigned i = 0; i < 10; ++i) {
    decisions.clear();

    // Heavy insert-delete workload
    for(unsigned nSteps = 0; nSteps < 1000; ++nSteps) {
      lastTreeGraph = tree.dumpGraphviz();

      // Decide whether to insert or remove a random item
      auto decisionFloat = TemplateMagic::random.getSingle<double>(0.0, 1.0);
      if(decisionFloat >= static_cast<double>(inTree.size()) / nKeys) {
        addElement(lastTreeGraph);
      } else {
        removeElement(lastTreeGraph);
      }

      fullValidation(lastTreeGraph);
    }

    auto matchingThroughIteration = TemplateMagic::zipMap(
      tree,
      inTree,
      [&](const unsigned& treeValue, const unsigned& testValue) -> bool {
        if(treeValue != testValue) {
          std::cout << "Expected " << testValue << ", got " << treeValue << std::endl;
          return false;
        }

        return true;
      }
    );

    BOOST_REQUIRE_MESSAGE(
      TemplateMagic::all_of(matchingThroughIteration),
      "BTree through-iteration does not yield the same elements as expected!\n"
        << tree.dumpGraphviz()
    );

    // Fill'er up all the way
    while(inTree.size() != nKeys) {
      std::string lastTreeGraph = tree.dumpGraphviz();

      addElement(lastTreeGraph);
      fullValidation(lastTreeGraph);
    }

    // Empty the tree
    while(!inTree.empty()) {
      std::string lastTreeGraph = tree.dumpGraphviz();

      removeElement(lastTreeGraph);
      fullValidation(lastTreeGraph);
    }
  }
}

BOOST_AUTO_TEST_CASE(setRemovalDefect) {
  /* Sets and maps cannot use std::remove_if! Not a defect. */

  auto testSet = TemplateMagic::moveIf(
    std::set<unsigned> {5, 2, 9, 1, 3, 4, 8},
    [](const auto& value) -> bool {
      return value % 2 != 0;
    }
  );

  BOOST_CHECK((testSet == std::set<unsigned> {1, 3, 5, 9}));
}

BOOST_AUTO_TEST_CASE(selectTestCases) {
  auto ragged2D = std::vector<
    std::vector<unsigned>
  > {
    {4, 9, 3},
    {2, 11, 6},
    {30, 4},
    {9, 44, 33, 12}
  };

  auto smallestVector = TemplateMagic::select(
    ragged2D,
    std::less<unsigned>(),
    [](const auto& group) -> unsigned {
      return group.size();
    }
  );

  BOOST_CHECK((*smallestVector == std::vector<unsigned> {30, 4}));

  auto largestVector = TemplateMagic::select(
    ragged2D,
    std::greater<unsigned>(),
    [](const auto& group) -> unsigned {
      return group.size();
    }
  );

  BOOST_CHECK((*largestVector == std::vector<unsigned> {9, 44, 33, 12}));

}

BOOST_AUTO_TEST_CASE(randomTests) {
  BOOST_CHECK(TemplateMagic::random.getSingle<unsigned>(0, 0) == 0);
}

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

struct Noisy {
  int f = -3;

  static unsigned nConstructed;

  Noisy() { ++nConstructed; }
  Noisy(const Noisy&) { ++nConstructed; }
  Noisy(Noisy&&) { ++nConstructed; }
};

unsigned Noisy::nConstructed = 0;

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

BOOST_AUTO_TEST_CASE(optionalEnhancementTests) {
  using namespace TemplateMagic;

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

  static_assert(
    std::is_same<
      decltype(fourth),
      detail::Optional<double>
    >::value,
    "No"
  );
}
