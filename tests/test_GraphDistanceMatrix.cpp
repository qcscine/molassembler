#include "BoostTestingHeader.h"
#include "GraphDistanceMatrix.h"

/*   Member listing
 * y 1 AdjacencyList constructor
 * y 2 getMatrixRef
 * y 3 altering operator ()
 * y 4 non-altering operator ()
 */

BOOST_AUTO_TEST_CASE( GraphDistanceMatrixTests ) {
  using namespace MoleculeManip;

  EdgeList edges {
    Edge(0, 1, BondType::Single),
    Edge(1, 2, BondType::Single),
    Edge(1, 4, BondType::Single),
    Edge(2, 3, BondType::Single),
    Edge(3, 4, BondType::Single),
    Edge(4, 5, BondType::Single),
    Edge(5, 6, BondType::Single),
    Edge(5, 7, BondType::Single)
  };

  /* AdjacencyMatrix(AdjacencyList(EdgeList))
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
  /* TODO change implementation to avoid "most vexing parse":
   * AdjacencyMatrix testInstance( ... );
   *
   * -> can be interpreted as function declaration! How to avoid?
   */
  auto testInstance = GraphDistanceMatrix(
    AdjacencyMatrix(
      AdjacencyList(
        edges
      )
    )
  );

  BOOST_CHECK(testInstance.N == 8);

  // original edges should be intact!
  for(const auto& edge: edges) {
    /* 4 */
    BOOST_CHECK(
      testInstance(
        edge.j, // since order shouldn't matter
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
