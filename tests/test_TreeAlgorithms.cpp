#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include "TreeAlgorithms.h"
#include <set>

using namespace MoleculeManip;
using namespace BasicTree;

template<typename T>
bool vectorsContainSameElements(
  const std::vector<T>& a,
  const std::vector<T>& b
) {
  std::set<T> sA(a.begin(), a.end()), sB(b.begin(), b.end());
  std::vector<T> difference;
  std::set_difference(
    sA.begin(),
    sA.end(),
    sB.begin(),
    sB.end(),
    std::back_inserter(difference)
  );
  return difference.size() == 0;
}

BOOST_AUTO_TEST_CASE( tree_output ) {
  std::shared_ptr<NodeType> rootPtr = std::make_shared<NodeType>(
    0
  );

  auto child = addChild<AtomIndexType>(rootPtr, 4);
  addChild<AtomIndexType>(child, 9);
  addChild<AtomIndexType>(child, 5);

  auto child2 = addChild<AtomIndexType>(rootPtr, 3);
  addChild<AtomIndexType>(child2, 6);

  std::cout << rootPtr << std::endl;

}

BOOST_AUTO_TEST_CASE( makeTreeTest ) {
  AdjacencyList test(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(0, 2, BondType::Single),
      Edge(1, 2, BondType::Single)
    })
  );

  /* should be expanded as
   * 0
   * |
   * 1
   * |
   * 2
   * |
   * 0
   *
   * 0-1-2-0
   *  `2
   */  

  auto treeStruct = makeTree(test);
  std::cout << treeStruct.rootPtr << std::endl;

}
