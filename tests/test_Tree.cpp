#include "BoostTestingHeader.h"
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

/* TODO
 * - it's obvious this API is shit... Having to notify the children of their
 *   parent shared_ptr ?? The fuck was I thinking
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
  auto secondChildPtr = firstInstancePtr -> addChild(9);
  secondChildPtr -> parentWeakPtr = firstInstancePtr;

  BOOST_CHECK(!secondChildPtr -> isRoot());
  BOOST_CHECK(secondChildPtr -> isLeaf());

  auto thirdChildPtr = firstInstancePtr -> addChild(11);
  thirdChildPtr -> parentWeakPtr = firstInstancePtr;

  BOOST_CHECK(!thirdChildPtr -> isRoot());
  BOOST_CHECK(thirdChildPtr -> isLeaf());
}
