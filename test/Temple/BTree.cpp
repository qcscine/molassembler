/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/constexpr/Array.h"
#include "Molassembler/Temple/constexpr/BTree.h"
#include "Molassembler/Temple/constexpr/Jsf.h"

#include <set>
#include <iostream>

using namespace Scine::Molassembler;

extern Temple::Generator<> generator;

namespace BTreeStaticTests {

constexpr Temple::BTree<unsigned, 3, 20> generateTree() {
  Temple::BTree<unsigned, 3, 20> tree;

  tree.insert(9);
  tree.insert(3);
  tree.insert(5);
  tree.insert(20);

  return tree;
}

constexpr auto testTree = generateTree();
static_assert(testTree.size() == 4, "Size of generated tree is wrong!");

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
  Temple::BTreeProperties::minHeight(5, 3) == 0
  && Temple::BTreeProperties::minHeight(35, 3) == 1
  && Temple::BTreeProperties::minHeight(215, 3) == 2,
  "minHeight function is wrong"
);

static_assert(
  Temple::BTreeProperties::maxNodesInTree(0, 3) == 1
  && Temple::BTreeProperties::maxNodesInTree(1, 3) == 7
  && Temple::BTreeProperties::maxNodesInTree(2, 3) == 43
  && Temple::BTreeProperties::maxNodesInTree(3, 3) == 259,
  "maxNodesInTree is wrong"
);

} // namespace BTreeStaticTests

inline unsigned popRandom(std::set<unsigned>& values) {
  auto it = values.begin();

  std::advance(
    it,
    Temple::Random::getSingle<unsigned>(0, values.size() - 1, generator.engine)
  );

  auto value = *it;

  values.erase(it);

  return value;
}

BOOST_AUTO_TEST_CASE(ConstexprBTree, *boost::unit_test::label("Temple")) {
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

  Temple::BTree<unsigned, 3, nKeys> tree;

  std::string lastTreeGraph;

  std::vector<std::string> decisions;

  auto addElement = [&](const std::string& treeGraph) {
    // Add an element
    auto toAdd = popRandom(notInTree);
    decisions.emplace_back("i"s + std::to_string(toAdd));

    BOOST_TEST_CONTEXT(
      "Element insertion threw. Operation sequence: "
      << Temple::condense(decisions)
      << ". Prior to last operation: \n"
      << treeGraph << "\n\n After last operation: \n"
      << tree.dumpGraphviz()
    ) {
      BOOST_REQUIRE_NO_THROW(tree.insert(toAdd));
    }

    inTree.insert(toAdd);
  };

  auto removeElement = [&](const std::string& treeGraph) {
    // Remove an element
    auto toRemove = popRandom(inTree);
    decisions.emplace_back("r"s + std::to_string(toRemove));

    BOOST_TEST_CONTEXT(
      "Tree element removal failed. Operation sequence: "
        << Temple::condense(decisions)
        << ". Prior to last operation: \n"
        << treeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    ) {
      BOOST_REQUIRE_NO_THROW(tree.remove(toRemove));
    }

    notInTree.insert(toRemove);
  };

  auto fullValidation = [&](const std::string& treeGraph) {
    // Validate the tree
    BOOST_TEST_CONTEXT(
      "Tree validation failed. Operation sequence: "
        << Temple::condense(decisions)
        << ". Prior to last operation: \n"
        << treeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    ) {
      BOOST_REQUIRE_NO_THROW(tree.validate());
    }

    // Check that elements that weren't inserted aren't falsely contained
    auto notInsertedButContained = Temple::copy_if(
      notInTree,
      [&](const auto& notInTreeValue) -> bool {
        return tree.contains(notInTreeValue);
      }
    );

    // Check that all elements are truly contained or not
    BOOST_REQUIRE_MESSAGE(
      notInsertedButContained.empty(),
      "Not all elements recorded as not in the tree are recognized as such!\n"
        << "Found in the tree, but should not be present: "
        << Temple::condense(notInsertedButContained)
        << "\nSequence of operations: "
        << Temple::condense(decisions)
        << ". Prior to last operation: \n"
        << treeGraph << "\n\n After last operation: \n"
        << tree.dumpGraphviz()
    );

    auto insertedNotContained = Temple::copy_if(
      inTree,
      [&](const auto& inTreeValue) -> bool {
        return !tree.contains(inTreeValue);
      }
    );

    BOOST_REQUIRE_MESSAGE(
      insertedNotContained.empty(),
      "Not all elements recorded as contained in the tree are recognized as such!\n"
        << "Not found in the tree: "
        << Temple::condense(insertedNotContained)
        << "\nSequence of operations: "
        << Temple::condense(decisions)
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
      auto decisionFloat = Temple::Random::getSingle<double>(0.0, 1.0, generator.engine);
      if(decisionFloat >= static_cast<double>(inTree.size()) / nKeys) {
        addElement(lastTreeGraph);
      } else {
        removeElement(lastTreeGraph);
      }

      fullValidation(lastTreeGraph);
    }

    BOOST_REQUIRE_MESSAGE(
      Temple::all_of(
        Temple::Adaptors::zip(tree, inTree),
        [&](const unsigned treeValue, const unsigned testValue) -> bool {
          if(treeValue != testValue) {
            std::cout << "Expected " << testValue << ", got " << treeValue << std::endl;
            return false;
          }

          return true;
        }
      ),
      "BTree through-iteration does not yield the same elements as expected!\n"
        << tree.dumpGraphviz()
    );

    // Fill'er up all the way
    while(inTree.size() != nKeys) {
      lastTreeGraph = tree.dumpGraphviz();

      addElement(lastTreeGraph);
      fullValidation(lastTreeGraph);
    }

    // Empty the tree
    while(!inTree.empty()) {
      lastTreeGraph = tree.dumpGraphviz();

      removeElement(lastTreeGraph);
      fullValidation(lastTreeGraph);
    }
  }
}

namespace {

/* Test that if a BTree is instantiated with a specific size, that size
 * definitely fits in the tree
 */

template<size_t minOrder, size_t nElements>
constexpr bool BTreeAllocatedSizeSufficient() {
  Temple::BTree<unsigned, minOrder, nElements> tree;

  for(unsigned i = 0; i < nElements; ++i) {
    tree.insert(i);
  }

  return true;
}

template<size_t minOrder, size_t ... nElements>
constexpr bool testAllBTrees(std::index_sequence<nElements...> /* elements */) {
  Temple::Array<bool, sizeof...(nElements)> results {{
    BTreeAllocatedSizeSufficient<minOrder, 5 + nElements>()...
  }};

  for(unsigned i = 0; i < sizeof...(nElements); ++i) {
    if(!results.at(i)) {
      return false;
    }
  }

  return true;
}

template<size_t ... minOrders>
constexpr bool testAllBTrees(std::index_sequence<minOrders...> /* elements */) {
  Temple::Array<bool, sizeof...(minOrders)> results {{
    testAllBTrees<2 + minOrders>(std::make_index_sequence<45>{})... // Test sizes 5->50
  }};

  for(unsigned i = 0; i < sizeof...(minOrders); ++i) {
    if(!results.at(i)) {
      return false;
    }
  }

  return true;
}

constexpr bool testAllBTrees() {
  return testAllBTrees(std::make_index_sequence<3>{}); // Test min orders 2->5
}

static_assert(
  testAllBTrees(),
  "For some B-Trees, you cannot fit as many elements in as requested at instantiation"
);

} // namespace
