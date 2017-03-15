#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TreeTests
#include <boost/test/unit_test.hpp>

#include <iostream>
#include "Tree.h"


/*       AdjacencyList member listing
 * t  #  ------------------------------
 * y  1  key constructor
 * y  2  addChild(Ptr)
 * y  3  addChild(key)
 * y  4  isRoot
 * y  5  isLeaf
 * y  6  operator ==
 */

/* TODO
 */

BOOST_AUTO_TEST_CASE( treeTests ) {
  using NodeType = Tree::Node<unsigned>;
  using namespace Tree;

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
  // tree is now 4 -> 5
  BOOST_CHECK(!firstInstancePtr -> isLeaf());

  // 4 (key addChild)
  auto secondChildPtr = firstInstancePtr -> addChild(9);
  // tree is now 4 -> 5 
  //               -> 9

  BOOST_CHECK(secondChildPtr->depth() == 1);

  BOOST_CHECK(!secondChildPtr -> isRoot());
  BOOST_CHECK(secondChildPtr -> isLeaf());

  auto thirdChildPtr = secondChildPtr -> addChild(11);
  /* tree is now 4 -> 5 
   *               -> 9 -> 11
   */

  BOOST_CHECK(thirdChildPtr->depth() == 2);

  BOOST_CHECK(!thirdChildPtr -> isRoot());
  BOOST_CHECK(thirdChildPtr -> isLeaf());
  BOOST_CHECK(!secondChildPtr -> isLeaf());

  // operator ==
  BOOST_CHECK(*firstInstancePtr != *secondChildPtr);
  BOOST_CHECK(*firstInstancePtr == *firstInstancePtr);

  // make an identical tree elsewhere
  auto identicalTreePtr = std::make_shared<NodeType>(4);
  auto levelOnePtr = identicalTreePtr -> addChild(9);
  levelOnePtr -> addChild(11);
  identicalTreePtr -> addChild(5);

  BOOST_CHECK(levelOnePtr->depth() == 1);

  BOOST_CHECK(*firstInstancePtr == *identicalTreePtr);

  // short-initialization
  auto shortInitializedTree = nodePtr(4u, {
    nodePtr(5u),
    nodePtr(9u, {11u})
  });

  BOOST_CHECK(*firstInstancePtr == *shortInitializedTree);
}
