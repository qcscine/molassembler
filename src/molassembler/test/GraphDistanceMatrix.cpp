#define BOOST_TEST_MODULE GraphDistanceMatrixTests
#include <boost/test/unit_test.hpp>

#include "Molecule.h"
#include "GraphDistanceMatrix.h"

// Testing help
#include "RepeatedElementCollection.h"

/*   Member listing
 * y 1 Molecule constructor
 * y 2 getMatrixRef
 * y 3 altering operator ()
 * y 4 non-altering operator ()
 */

BOOST_AUTO_TEST_CASE( GraphDistanceMatrixTests ) {
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

  /* AdjacencyMatrix(Molecule(Edge))
   *
   *    0 1 2 3 4 5 6 7
   *  0 · 1 0 0 0 0 0 0
   *  1 · · 1 0 1 0 0 0
   *  2 · · · 1 0 0 0 0
   *  3 · · · · 1 0 0 0
   *  4 · · · · · 1 0 0
   *  5 · · · · · · 1 1
   *  6 · · · · · · · 0
   *  7 · · · · · · · ·
   *
   * GraphDistanceMatrix(Ans)
   *
   *    0 1 2 3 4 5 6 7
   *  0 · 1 2 3 2 3 4 4
   *  1 · · 1 2 1 2 3 3
   *  2 · · · 1 2 3 4 4
   *  3 · · · · 1 2 3 3
   *  4 · · · · · 1 2 2
   *  5 · · · · · · 1 1
   *  6 · · · · · · · 2
   *  7 · · · · · · · ·
   *
   */

  /* 1 */
  auto testInstance = GraphDistanceMatrix {
    AdjacencyMatrix {
      Molecule { std::move(graph) }
    }
  };

  BOOST_CHECK(testInstance.N == 8);

  // original edges should be intact!
  for(const auto& edge: edges) {
    /* 4 */
    BOOST_CHECK(
      testInstance(
        edge.j,
        edge.i
      ) == 1
    );
  }

  // check expectations for the hardcoded example from above
  auto checkRow = [&testInstance](
    const unsigned& rowNumber,
    const std::vector<unsigned>& expectedValues
  ) {
    for(unsigned i = 0; i < expectedValues.size(); i++) {
      BOOST_CHECK(
        testInstance(
          rowNumber,
          rowNumber + i + 1
        ) == expectedValues[i]
      );
    }
  };

  checkRow(0, {1, 2, 3, 2, 3, 4, 4});
  checkRow(1, {1, 2, 1, 2, 3, 3});
  checkRow(2, {1, 2, 3, 4, 4});
  checkRow(3, {1, 2, 3, 3});
  checkRow(4, {1, 2, 2});
  checkRow(5, {1, 1});
  checkRow(6, {2});

}
