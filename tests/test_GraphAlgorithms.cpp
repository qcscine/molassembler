#include "BoostTestingHeader.h"

#include "AdjacencyListAlgorithms.h"
#include "StdlibTypeAlgorithms.h"


/* Connected components tests */
BOOST_AUTO_TEST_CASE( GraphAlgorithmsTests ) {
  using namespace MoleculeManip;

  { 
    AdjacencyList test(
      EdgeList({
        Edge(0, 1, BondType::Single),
        Edge(1, 2, BondType::Single),
        Edge(2, 3, BondType::Single)
      })
    );

    BOOST_CHECK(AdjacencyListAlgorithms::numConnectedComponents(test) == 1);
  }

  {
    AdjacencyList test(
      EdgeList({
        Edge(0, 1, BondType::Single),
        Edge(0, 2, BondType::Single),
        Edge(1, 2, BondType::Single)
      })
    );

    BOOST_CHECK(AdjacencyListAlgorithms::numConnectedComponents(test) == 1);
  }

  {
    AdjacencyList test(
      EdgeList({
        Edge(0, 1, BondType::Single),
        Edge(1, 2, BondType::Single),
        Edge(3, 4, BondType::Single)
      })
    );

    BOOST_CHECK(AdjacencyListAlgorithms::numConnectedComponents(test) == 2);
  }

  /* Cycle detection tests */
  {
    AdjacencyList test(
      EdgeList({
        Edge(0, 1, BondType::Single),
        Edge(0, 2, BondType::Single),
        Edge(1, 2, BondType::Single)
      })
    );

    auto result = AdjacencyListAlgorithms::detectCycles(test);

    BOOST_CHECK(result.size() == 1);
    BOOST_CHECK(
      StdlibTypeAlgorithms::vectorToSet(result[0]) 
      == std::set<AtomIndexType>({0, 1, 2})
    );
  }

}
