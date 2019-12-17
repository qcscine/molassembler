/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "temple/BoundedNodeTrie.h"
#include "temple/OrderedPair.h"
#include "temple/Poset.h"
#include "temple/Random.h"
#include "temple/Stringify.h"
#include "temple/constexpr/Jsf.h"

#include <algorithm>
#include <random>

extern temple::jsf::Generator<> generator;

BOOST_AUTO_TEST_CASE(OrderedPairTests) {
  temple::OrderedPair<unsigned> a {14u, 3u};
  BOOST_CHECK(a.front() < a.back());
  BOOST_CHECK(a.first < a.second);
  BOOST_CHECK(
    std::is_sorted(std::begin(a), std::end(a))
  );
}

BOOST_AUTO_TEST_CASE(PosetTests) {
  constexpr unsigned nTests = 10;
  constexpr unsigned nValues = 20;
  constexpr unsigned maxValue = 40;

  auto compareFirstDigit = [](const unsigned a, const unsigned b) -> bool {
    return (a / 10) < (b / 10);
  };

  std::vector<unsigned> unorderedValues;
  for(unsigned testNum = 0; testNum < nTests; ++testNum) {
    unorderedValues = temple::random::getN<unsigned>(0, maxValue, nValues, generator.engine);
    temple::Poset<unsigned> f {unorderedValues};

    f.orderUnordered(compareFirstDigit);
    f.orderUnordered(std::less<>());
    f.finalize();

    BOOST_CHECK(
      std::is_sorted(
        std::begin(f),
        std::end(f),
        [](const auto& a, const auto& b) -> bool {
          return std::less<>()(
            a.values.front(),
            b.values.front()
          );
        }
      )
    );
  }
}

template<
  typename ChoiceIndex,
  typename Generator
>
struct ChooseFunctor {
  Generator& prngRef;

  ChooseFunctor(Generator& prng) : prngRef(prng) {}

  int euclideanModulo(const int a, const int base) {
    return ((a % base) + base) % base;
  }

  int distance(const int i, const int j, const int U) {
    assert(i >= 0 && j >= 0 && U > 1);

    return std::min(
      euclideanModulo(i - j, U),
      euclideanModulo(j - i, U)
    );
  }

  double merit(const ChoiceIndex choice, const boost::dynamic_bitset<>& childExists) {
    double value = 0.0;

    const ChoiceIndex U = childExists.size();

    for(ChoiceIndex other = 0; other < U; ++other) {
      if(other == choice) {
        continue;
      }

      if(childExists.test(other)) {
        int choiceDistance = distance(choice, other, U);
        assert(choiceDistance >= 0);
        assert(choiceDistance < std::numeric_limits<ChoiceIndex>::max());
        value += choiceDistance;
      }
    }

    assert(value >= 0);

    return value;
  }

  ChoiceIndex operator () (const std::vector<ChoiceIndex>& viableChildren, const boost::dynamic_bitset<>& existingChildren) {
    assert(!viableChildren.empty());

    /* In case not all children already exist, it makes sense to calculate
     * merits for all choices at this level.
     */
    if(!existingChildren.all()) {
      std::vector<ChoiceIndex> bestChoices;
      double bestMerit = 0;
      for(ChoiceIndex viableChild : viableChildren) {
        if(!existingChildren.test(viableChild)) {
          double choiceMerit = merit(viableChild, existingChildren);
          if(choiceMerit > bestMerit) {
            bestChoices = {viableChild};
            bestMerit = choiceMerit;
          } else if(choiceMerit == bestMerit) {
            bestChoices.push_back(viableChild);
          }
        }
      }

      assert(!bestChoices.empty());

      return bestChoices.at(
        std::uniform_int_distribution<ChoiceIndex>(
          0,
          bestChoices.size() - 1
        )(prngRef.engine)
      );
    }

    /* All children exist, so no point in calculating merits: we choose from
     * the viable ones
     */
    return viableChildren.at(
      std::uniform_int_distribution<ChoiceIndex>(0, viableChildren.size() - 1)(prngRef.engine)
    );
  }
};

template<typename ChoiceIndex, typename Generator>
auto make_ChooseFunctor(
  const temple::BoundedNodeTrie<ChoiceIndex>& /* trie */,
  Generator& prng
) {
  return ChooseFunctor<ChoiceIndex, Generator>(prng);
}


BOOST_AUTO_TEST_CASE(asdf) {
  {
    std::vector<std::uint8_t> boundaries {4, 2, 3};

    using TrieType = temple::BoundedNodeTrie<std::uint8_t>;

    TrieType trie {boundaries};

    std::vector<TrieType::ChoiceList> sampleValueLists {
      {0, 0, 0},
      {0, 1, 2},
      {3, 0, 1}
    };

    BOOST_CHECK(trie.capacity() == 24);

    for(auto& sampleValueList : sampleValueLists) {
      BOOST_CHECK(!trie.contains(sampleValueList));
    }

    for(unsigned i = 0; i < sampleValueLists.size(); ++i) {
      trie.insert(sampleValueLists.at(i));

      for(unsigned j = 0; j <= i; ++j) {
        BOOST_CHECK(trie.contains(sampleValueLists.at(j)));
      }
      for(unsigned j = i + 1; j < sampleValueLists.size(); ++j) {
        BOOST_CHECK(!trie.contains(sampleValueLists.at(j)));
      }

      BOOST_CHECK(trie.size() == i + 1);
    }
  }
  {
    using TrieType = temple::BoundedNodeTrie<std::uint8_t>;
    TrieType::ChoiceList boundaries;
    boundaries.push_back(4);

    TrieType trie {boundaries};

    std::vector<TrieType::ChoiceList> sampleValueLists {
      {0},
      {3},
      {2}
    };

    BOOST_CHECK(trie.capacity() == 4);

    for(auto& sampleValueList : sampleValueLists) {
      BOOST_CHECK(!trie.contains(sampleValueList));
    }

    for(unsigned i = 0; i < sampleValueLists.size(); ++i) {
      trie.insert(sampleValueLists.at(i));

      for(unsigned j = 0; j <= i; ++j) {
        BOOST_CHECK(trie.contains(sampleValueLists.at(j)));
      }
      for(unsigned j = i + 1; j < sampleValueLists.size(); ++j) {
        BOOST_CHECK(!trie.contains(sampleValueLists.at(j)));
      }

      BOOST_CHECK(trie.size() == i + 1);
    }
  }
  {
    // Generate new list every time test
    using TrieType = temple::BoundedNodeTrie<std::uint8_t>;
    TrieType::ChoiceList boundaries {6, 4, 2};
    TrieType trie {boundaries};
    auto chooseFunctor = make_ChooseFunctor(trie, generator);

    std::set<TrieType::ChoiceList> listsSet;

    TrieType::ChoiceList list;
    while(trie.size() != trie.capacity()) {
      list = trie.generateNewEntry(chooseFunctor);

      BOOST_REQUIRE(listsSet.count(list) == 0);
      listsSet.insert(list);
    }
  }
}
