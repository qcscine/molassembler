#include <boost/test/results_collector.hpp>
#include <boost/test/unit_test.hpp>

#include "temple/BTree.h"
#include "temple/Random.h"
#include "temple/Containers.h"

#include <set>
#include <iostream>

extern temple::Generator prng;

inline bool lastTestPassed() {
  using namespace boost::unit_test;

  test_case::id_t id = framework::current_test_case().p_id;
  test_results rez = results_collector.results(id);
  return rez.passed();
}

inline unsigned popRandom(std::set<unsigned>& values) {
  auto it = values.begin();

  std::advance(
    it,
    prng.getSingle<unsigned>(0, values.size() - 1)
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

  temple::nonconstexpr::BTree<unsigned, 3> tree;

  std::string lastTreeGraph;

  std::vector<std::string> decisions;

  auto addElement = [&](const std::string& treeGraph) {
    // Add an element
    auto toAdd = popRandom(notInTree);
    decisions.emplace_back("i"s + std::to_string(toAdd));

    BOOST_CHECK_NO_THROW(tree.insert(toAdd));
    BOOST_REQUIRE_MESSAGE(
      lastTestPassed(),
      "Element insertion failed. Operation sequence: "
        << temple::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << treeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    inTree.insert(toAdd);
  };

  auto removeElement = [&](const std::string& treeGraph) {
    // Remove an element
    auto toRemove = popRandom(inTree);
    decisions.emplace_back("r"s + std::to_string(toRemove));

    BOOST_CHECK_NO_THROW(tree.remove(toRemove));
    BOOST_REQUIRE_MESSAGE(
      lastTestPassed(),
      "Tree element removal failed. Operation sequence: "
        << temple::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << treeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    notInTree.insert(toRemove);
  };

  auto fullValidation = [&](const std::string& treeGraph) {
    // Validate the tree
    BOOST_CHECK_NO_THROW(tree.validate());
    BOOST_REQUIRE_MESSAGE(
      lastTestPassed(),
      "Tree validation failed. Operation sequence: "
        << temple::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << treeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    auto notInsertedNotContained = temple::mapToVector(
      notInTree,
      [&](const auto& notInTreeValue) -> bool {
        return !tree.contains(notInTreeValue);
      }
    );

    // Check that all elements are truly contained or not
    BOOST_REQUIRE_MESSAGE(
      temple::all_of(notInsertedNotContained),
      "Not all elements recorded as not in the tree are recognized as such!\n"
        << "Found in the tree, but should not be present: "
        << temple::condenseIterable(
          temple::moveIf(
            temple::zipMap(
              notInsertedNotContained,
              notInTree,
              [](const bool passed, const unsigned value) -> std::string {
                if(!passed) {
                  return std::to_string(value);
                }

                return "";
              }
            ),
            [](const std::string& str) -> bool {
              return !str.empty();
            }
          )
        ) << "\nSequence of operations: "
        << temple::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << treeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    auto insertedContained = temple::mapToVector(
      inTree,
      [&](const auto& inTreeValue) -> bool {
        return tree.contains(inTreeValue);
      }
    );

    BOOST_REQUIRE_MESSAGE(
      temple::all_of(insertedContained),
      "Not all elements recorded as contained in the tree are recognized as such!\n"
        << "Not found in the tree: "
        << temple::condenseIterable(
          temple::moveIf(
            temple::zipMap(
              insertedContained,
              inTree,
              [](const bool passed, const unsigned value) -> std::string {
                if(!passed) {
                  return std::to_string(value);
                }

                return "";
              }
            ),
            [](const std::string& str) -> bool {
              return !str.empty();
            }
          )
        ) << "\nSequence of operations: "
        << temple::condenseIterable(decisions)
        << ". Prior to last operation: \n"
        << treeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );
  };

  for(unsigned i = 0; i < 10; ++i) {
    decisions.clear();

    // Heavy insert-delete workload
    for(unsigned nSteps = 0; nSteps < 1000; ++nSteps) {
      lastTreeGraph = tree.dumpGraphviz();

      // Decide whether to insert or remove a random item
      auto decisionFloat = prng.getSingle<double>(0.0, 1.0);
      if(decisionFloat >= static_cast<double>(inTree.size()) / nKeys) {
        addElement(lastTreeGraph);
      } else {
        removeElement(lastTreeGraph);
      }

      fullValidation(lastTreeGraph);
    }

    auto matchingThroughIteration = temple::zipMap(
      tree,
      inTree,
      [&](const unsigned treeValue, const unsigned testValue) -> bool {
        if(treeValue != testValue) {
          std::cout << "Expected " << testValue << ", got " << treeValue << std::endl;
          return false;
        }

        return true;
      }
    );

    BOOST_REQUIRE_MESSAGE(
      temple::all_of(matchingThroughIteration),
      "BTree through-iteration does not yield the same elements as expected!\n"
        << tree.dumpGraphviz()
    );

    // Fill'er up all the way
    while(inTree.size() != nKeys) {
      std::string anotherTreeGraph = tree.dumpGraphviz();

      addElement(anotherTreeGraph);
      fullValidation(anotherTreeGraph);
    }

    // Empty the tree
    while(!inTree.empty()) {
      std::string anotherTreeGraph = tree.dumpGraphviz();

      removeElement(anotherTreeGraph);
      fullValidation(anotherTreeGraph);
    }
  }
}
