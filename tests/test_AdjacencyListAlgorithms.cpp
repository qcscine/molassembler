#include "BoostTestingHeader.h"

#include "AdjacencyListAlgorithms.h"
#include "StdlibTypeAlgorithms.h"

/* TODO
 * - A potential problem: tree generation can change according to the order 
 *   that nodes of the AdjacencyList graph are visited in by BFS. To test, use
 *   a five ring variant of the existing adjacencylist and make random
 *   permutations of the order of adjacencylist addition (without Edges 
 *   intermediate, it orders it's parameters)
 */


/*       AdjacencyListAlgorithms listing
 * t  #  -------------------------------
 * y  1  _connectedComponents
 * y  2  numConnectedComponents
 * y  3  connectedComponentGroups
 *    4  detectCycles
 */

BOOST_AUTO_TEST_CASE( adjacencyListAlgorithms ) {
  using namespace MoleculeManip;
  using namespace AdjacencyListAlgorithms;

  auto testInstance = AdjacencyList(
    Edges({
      {{0, 1}, BondType::Single},
      {{1, 2}, BondType::Single},
      {{1, 4}, BondType::Single},
      {{2, 3}, BondType::Single},
      {{3, 4}, BondType::Single},
      {{4, 5}, BondType::Single},
      {{5, 6}, BondType::Single},
      {{5, 7}, BondType::Single}
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
    StdlibTypeAlgorithms::vectorToSet(groupsVectors) 
    == StdlibTypeAlgorithms::vectorToSet(expectedGroups)
  );

  std::vector<AtomIndexType> BFSVisitSequence;
  BFSVisit(
    testInstance,
    0,
    [&BFSVisitSequence](const AtomIndexType& index) -> bool {
      BFSVisitSequence.push_back(index);
      return true;
    }
  );
  BOOST_CHECK(
    BFSVisitSequence 
    == std::vector<AtomIndexType>({
      0, 1, 2, 4, 3, 5, 6, 7
    })
  );

  std::vector<AtomIndexType> DFSVisitSequence;
  DFSVisit(
    testInstance,
    0,
    [&DFSVisitSequence](const AtomIndexType& index) -> bool {
      DFSVisitSequence.push_back(index);
      return true;
    }
  );
  BOOST_CHECK(
    DFSVisitSequence 
    == std::vector<AtomIndexType>({
      0, 1, 4, 5, 7, 6, 3, 2
    })
  );

  { // expansion
    AdjacencyList test(
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 2}, BondType::Single}
      })
    );

    auto treePtr = makeTree(test);

    auto comparisonTreePtr = Tree::nodePtr(0u, {
      Tree::nodePtr(1u, {2u}),
      Tree::nodePtr(2u)
    });

    BOOST_CHECK(*treePtr == *comparisonTreePtr);
  }
  { // BFS, DFS testing
    AdjacencyList test(
      Edges({
        {{0, 1}, BondType::Single},
        {{1, 2}, BondType::Single},
        {{1, 4}, BondType::Single},
        {{2, 3}, BondType::Single},
        {{3, 4}, BondType::Single},
        {{4, 5}, BondType::Single},
        {{5, 6}, BondType::Single},
        {{5, 7}, BondType::Single}
      })
    );

    auto treePtr = makeTree(test, 0);

    auto comparisonTreePtr = Tree::nodePtr(0u, {
      Tree::nodePtr(1u, {
        Tree::nodePtr(2u, {3u}),
        Tree::nodePtr(4u, {
          Tree::nodePtr(3u),
          Tree::nodePtr(5u, {6u, 7u})
        })
      })
    });

    BOOST_CHECK(*treePtr == *comparisonTreePtr);
  }
  { // order independence of tree generation on sequence in AdjacencyList

  }
}
