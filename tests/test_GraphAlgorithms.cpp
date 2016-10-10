#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include "GraphAlgorithms.h"

using namespace MoleculeManip;

BOOST_AUTO_TEST_CASE( chain ) {
  AdjacencyList test(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(1, 2, BondType::Single),
      Edge(2, 3, BondType::Single)
    })
  );

  BOOST_CHECK(GraphAlgorithms::num_connected_components(test) == 1);
}

BOOST_AUTO_TEST_CASE( cycle ) {
  AdjacencyList test(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(0, 2, BondType::Single),
      Edge(1, 2, BondType::Single)
    })
  );

  BOOST_CHECK(GraphAlgorithms::num_connected_components(test) == 1);
}

BOOST_AUTO_TEST_CASE( simple_disconnect ) {
  AdjacencyList test(
    EdgeList({
      Edge(0, 1, BondType::Single),
      Edge(1, 2, BondType::Single),
      Edge(3, 4, BondType::Single)
    })
  );

  BOOST_CHECK(GraphAlgorithms::num_connected_components(test) == 2);
}
