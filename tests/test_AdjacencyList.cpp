#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "AdjacencyList.h"

using namespace MoleculeManip;

/*       AdjacencyList member listing
 * t  #  ------------------------------
 *    1  Default constructor
 * y  2  Constructor from EdgeList
 * y  3  addSlot()
 * y  4  addAdjacency()
 *    5  clear()
 * y  6  removeAdjacency()
 *    7  getAdjacencies()
 * y  8  isAdjacent()
 * y  9  size()
 * y  10 validate()
 *    11 operator[]
 */

BOOST_AUTO_TEST_CASE(init_from_EdgeList) {
  EdgeList edgeList({
    Edge(0, 1, BondType::Single),
    Edge(1, 2, BondType::Single),
    Edge(1, 4, BondType::Single),
    Edge(2, 3, BondType::Single),
    Edge(3, 4, BondType::Single),
    Edge(4, 5, BondType::Single),
    Edge(5, 6, BondType::Single),
    Edge(5, 7, BondType::Single)
  });

  /* 2 */
  AdjacencyList testList(
    edgeList
  );

  /* 10 */
  BOOST_REQUIRE(testList.validate());

  /* 9 */
  BOOST_CHECK(testList.size() == 7); // is maximum AtomIndexType supplied here

  for(const auto& edge: edgeList) {
    /* 8 */
    BOOST_CHECK(testList.isAdjacent(edge.i, edge.j));
  }

  /* 4 */
  testList.addAdjacency(7, 0);
  BOOST_CHECK(testList.size() == 7); // size is unchanged

  /* 3 */
  testList.addSlot();
  BOOST_CHECK(testList.size() == 8); // size changes
  BOOST_REQUIRE(testList.validate());

  /* 6 */
  testList.removeAdjacency(2, 3);
  BOOST_CHECK(testList.size() == 8); // size unchanged
  BOOST_REQUIRE(testList.validate());

}
