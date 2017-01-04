#include "BoostTestingHeader.h"
#include "EdgeList.h"

/* List of EdgeList functions
 *
 * y  1  default constructor
 * y  2  constructor from vector of edges
 * y  3  add
 * y  4  clear
 * y  5  remove (edge)
 * y  6  remove (a, b)
 * y  7  get
 * y  8  isOrdered
 * y  9  size
 * y  10 begin-end iteration
 */

BOOST_AUTO_TEST_CASE( edgelist_all ) {
  using namespace MoleculeManip;

  std::vector<Edge> edges {
    Edge(0, 1, BondType::Single),
    Edge(1, 2, BondType::Single),
    Edge(1, 4, BondType::Single),
    Edge(2, 3, BondType::Single),
    Edge(3, 4, BondType::Single),
    Edge(4, 5, BondType::Single),
    Edge(5, 6, BondType::Single),
    Edge(5, 7, BondType::Single)
  };

  /* 2 */
  EdgeList testInstance(edges);

  /* 9 */
  BOOST_CHECK(testInstance.size() == 8);
  /* 8 */
  BOOST_REQUIRE(testInstance.isOrdered());

  /* 10 */
  // Ensure all edges in testInstance are also in the original vector of edges
  BOOST_CHECK(
    std::accumulate(
      testInstance.begin(),
      testInstance.end(),
      true,
      [&edges](const bool& carry, const auto& edge) {
        return (
          carry
          && std::find(
            edges.begin(),
            edges.end(),
            edge
          ) != edges.end()
        );
      }
    )
  );

  // And reverse
  BOOST_CHECK(
    std::accumulate(
      edges.begin(),
      edges.end(),
      true,
      [&testInstance](const bool& carry, const auto& edge) {
        return (
          carry
          && std::find(
            testInstance.begin(),
            testInstance.end(),
            edge
          ) != testInstance.end()
        );
      }
    )
  );

  /* 7 */
  BOOST_CHECK(
    testInstance.get(9, 0)
    == boost::optional<Edge>()
  );
  BOOST_CHECK(
    testInstance.get(5, 4)
    != boost::optional<Edge>()
  );


  /* 4 */
  testInstance.clear();
  BOOST_REQUIRE(testInstance.size() == 0);

  /* 3 */
  testInstance.add(
    Edge(0, 5, BondType::Single)
  );
  testInstance.add(
    Edge(3, 4, BondType::Single)
  );
  BOOST_REQUIRE(testInstance.size() == 2);

  /* 5, 6 */
  testInstance.remove(3, 4);
  testInstance.remove(
    Edge(0, 5, BondType::Single)
  );
  BOOST_REQUIRE(testInstance.size() == 0);

  /* 1 */
  // Overwrite with default constructor
  testInstance = EdgeList();

}
