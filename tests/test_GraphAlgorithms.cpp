#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>


#include "GraphAlgorithms.h"

using namespace MoleculeManip;

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

/* Connected components tests */
BOOST_AUTO_TEST_CASE( connected_components_chain ) {
  AdjacencyList test(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(1, 2, BondType::Single),
      Edge(2, 3, BondType::Single)
    })
  );

  BOOST_CHECK(GraphAlgorithms::numConnectedComponents(test) == 1);
}

BOOST_AUTO_TEST_CASE( connected_components_cycle ) {
  AdjacencyList test(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(0, 2, BondType::Single),
      Edge(1, 2, BondType::Single)
    })
  );

  BOOST_CHECK(GraphAlgorithms::numConnectedComponents(test) == 1);
}

BOOST_AUTO_TEST_CASE( connected_components_simple_disconnect ) {
  AdjacencyList test(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(1, 2, BondType::Single),
      Edge(3, 4, BondType::Single)
    })
  );

  BOOST_CHECK(GraphAlgorithms::numConnectedComponents(test) == 2);
}

/* Cycle detection tests */
BOOST_AUTO_TEST_CASE( cycle_detection_triangle ) {
  AdjacencyList test(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(0, 2, BondType::Single),
      Edge(1, 2, BondType::Single)
    })
  );

  /* is expanded as:
   * 0-1-2-0
   */  

  auto result = GraphAlgorithms::detectCycles(test);
  BOOST_CHECK(result.size() == 1);
  BOOST_CHECK(vectorsContainSameElements(
    result[0],
    std::vector<AtomIndexType>({0, 1, 2})
  ));
}
