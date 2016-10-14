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

  std::shared_ptr<NodeType> child = std::make_shared<NodeType>(
    rootPtr,
    1
  );

  std::shared_ptr<NodeType> child2 = std::make_shared<NodeType>(
    rootPtr,
    2
  );

  rootPtr -> addChild(child);
  rootPtr -> addChild(child2);

  auto childChild = std::make_shared<NodeType>(
    rootPtr,
    3
  );

  child -> addChild(childChild);


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
   */  

  auto treeStruct = makeTree(test);
  std::cout << treeStruct.rootPtr << std::endl;

}
