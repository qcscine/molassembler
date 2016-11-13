#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "AdjacencyListAlgorithms.h"

using namespace MoleculeManip;

/*       AdjacencyListAlgorithms listing
 * t  #  -------------------------------
 *    1  _connectedComponents
 *    2  numConnectedComponents
 *    3  connectedComponentGroups
 *    4  detectCycles
 */

// MORE AdjacencyList instances to test with
//
using namespace AdjacencyListAlgorithms;

template<typename T>
bool vectorToSetEquals(
  const std::vector<T>& a,
  const std::vector<T>& b
) {
  std::set<T> setA(a.begin(), a.end());
  std::set<T> setB(b.begin(), b.end());

  return std::equal(
    setA.begin(),
    setA.end(),
    setB.begin(),
    setB.end()
  );
}

BOOST_AUTO_TEST_CASE( adjacencyListAlgorithms ) {
  auto testInstance = AdjacencyList(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(1, 2, BondType::Single),
      Edge(1, 4, BondType::Single),
      Edge(2, 3, BondType::Single),
      Edge(3, 4, BondType::Single),
      Edge(4, 5, BondType::Single),
      Edge(5, 6, BondType::Single),
      Edge(5, 7, BondType::Single)
    })
  );

  auto countVector = _connectedComponents(testInstance);
  std::vector<unsigned> expected {1, 1, 1, 1, 1, 1, 1, 1};

  BOOST_CHECK(
    std::equal(
      countVector.begin(),
      countVector.end(),
      expected.begin(),
      expected.end()
    )
  );


  BOOST_CHECK(
    numConnectedComponents(testInstance) == 1
  );

  auto groupsVectors = connectedComponentGroups(testInstance);
  std::vector<
    std::vector<AtomIndexType>
  > expectedGroups {
    {0, 1, 2, 3, 4, 5, 6, 7}
  };

  BOOST_CHECK(
    vectorToSetEquals(
      groupsVectors,
      expectedGroups
    )
  );

}
