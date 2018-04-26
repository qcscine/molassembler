#define BOOST_TEST_MODULE AdjacencyMatrixTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "AdjacencyMatrix.h"

// Helper header
#include "RepeatedElementCollection.h"

/*   Member listing
 * y 1 Molecule constructor
 * y 2 getMatrixRef
 * y 3 altering operator ()
 * y 4 non-altering operator ()
 */

BOOST_AUTO_TEST_CASE( AdjacencyMatrix_all ) {
  using namespace molassembler;

  struct Edge {
    AtomIndexType i, j;
    BondType bty;

    Edge() = default;
    Edge(AtomIndexType i, AtomIndexType j, BondType bty)
      : i {i}, j {j}, bty {bty}
    {}
  };

  unsigned size = 8;
  std::vector<Edge> edges {
    {0, 1, BondType::Single},
    {1, 2, BondType::Single},
    {1, 4, BondType::Single},
    {2, 3, BondType::Single},
    {3, 4, BondType::Single},
    {4, 5, BondType::Single},
    {5, 6, BondType::Single},
    {5, 7, BondType::Single}
  };

  GraphType graph(size);

  for(const auto& edge : edges) {
    auto addPair = boost::add_edge(edge.i, edge.j, graph);
    graph[addPair.first].bondType = edge.bty;
  }

  for(AtomIndexType i = 0; i < size; ++i) {
    graph[i].elementType = Delib::ElementType::H;
  }

  /* 1 */
  // use uniform initialization syntax to avoid most vexing parse
  AdjacencyMatrix testInstance {
    Molecule { std::move(graph) }
  };

  BOOST_CHECK(testInstance.N == 8);

  for(const auto& edge: edges) {
    /* 4 */
    BOOST_CHECK(
      testInstance(
        edge.j,
        edge.i
      ) // the positions are boolean already
    );
  }

  /* 3 */
  testInstance(5, 2) = true;

  /* 2 */
  BOOST_CHECK(testInstance.getMatrixRef()(2, 5));
  // This is faulty, we cannot say anything about the state of the lower matrix
  // BOOST_CHECK(testInstance.getMatrixRef()(5, 2) == 0);
}
