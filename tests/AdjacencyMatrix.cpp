#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AdjacencyMatrixTests
#include <boost/test/unit_test.hpp>

#include "AdjacencyMatrix.h"

/*   Member listing
 * y 1 AdjacencyList constructor
 * y 2 getMatrixRef
 * y 3 altering operator ()
 * y 4 non-altering operator ()
 */

BOOST_AUTO_TEST_CASE( AdjacencyMatrix_all ) {
  using namespace MoleculeManip;
  Edges edges{
    {{0, 1}, BondType::Single},
    {{1, 2}, BondType::Single},
    {{1, 4}, BondType::Single},
    {{2, 3}, BondType::Single},
    {{3, 4}, BondType::Single},
    {{4, 5}, BondType::Single},
    {{5, 6}, BondType::Single},
    {{5, 7}, BondType::Single}
  };

  /* 1 */
  // use uniform initialization syntax to avoid most vexing parse
  AdjacencyMatrix testInstance {
    AdjacencyList {
      edges
    }
  };

  BOOST_CHECK(testInstance.N == 8);

  for(const auto& edge: edges) {
    /* 4 */
    BOOST_CHECK(
      testInstance(
        edge.first.second, // since order shouldn't matter
        edge.first.first
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
