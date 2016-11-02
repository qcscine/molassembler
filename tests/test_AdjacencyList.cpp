#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "AdjacencyList.h"

using namespace MoleculeManip;

BOOST_AUTO_TEST_CASE(asdf) {
  AdjacencyList testList(
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

  std::cout << testList.distancesMatrix() << std::endl;
}
