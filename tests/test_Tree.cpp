#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "Tree.h"


/*       AdjacencyList member listing
 * t  #  ------------------------------
 *    1  key constructor
 *    2  addChild(Ptr)
 *    3  addChild(key)
 *    4  isRoot
 *    5  isLeaf
 */

BOOST_AUTO_TEST_CASE( treeTests ) {
  using NodeType = BasicTree::Node<unsigned>;

  // 1
  auto firstInstancePtr = std::make_shared<NodeType>(4);
  BOOST_CHECK(firstInstancePtr -> key == 4);
  // 5, 6
  BOOST_CHECK(firstInstancePtr -> isRoot());
  BOOST_CHECK(firstInstancePtr -> isLeaf());

  // 3 (member addChild)
  firstInstancePtr -> addChild(
    std::make_shared<NodeType>(5)
  );
  BOOST_CHECK(!firstInstancePtr -> isLeaf());

  // 4 (key addChild)
  auto secondChild = firstInstancePtr -> addChild(9);
  BOOST_CHECK(!secondChild -> isRoot());
  BOOST_CHECK(secondChild -> isLeaf());

  auto thirdChild = firstInstancePtr -> addChild(11);
  BOOST_CHECK(!thirdChild -> isRoot());
  BOOST_CHECK(thirdChild -> isLeaf());
}
