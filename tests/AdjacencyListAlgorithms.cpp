#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AdjacencyListAlgorithmsTests
#include <boost/test/unit_test.hpp>

#include "AdjacencyListAlgorithms.h"
#include "StdlibTypeAlgorithms.h"

// Helper header
#include "RepeatedElementCollection.h"

/* TODO
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
    makeRepeatedElementCollection(Delib::ElementType::H, 8),
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
  { // test BFS / DFS without depth limitation
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
  }
  { // test BFS / DFS with depth limitation
    std::vector<AtomIndexType> BFSVisitSequence;
    BFSVisit(
      testInstance,
      0,
      [&BFSVisitSequence](const AtomIndexType& index) -> bool {
        BFSVisitSequence.push_back(index);
        return true;
      },
      4
    );
    BOOST_CHECK(
      BFSVisitSequence 
      == std::vector<AtomIndexType>({
        0, 1, 2, 4, 3, 5 
      })
    );

    std::vector<AtomIndexType> DFSVisitSequence;
    DFSVisit(
      testInstance,
      0,
      [&DFSVisitSequence](const AtomIndexType& index) -> bool {
        DFSVisitSequence.push_back(index);
        return true;
      },
      4
    );
    BOOST_CHECK(
      DFSVisitSequence 
      == std::vector<AtomIndexType>({
        0, 1, 4, 5, 3, 2
      })
    );
  }
  auto testExpansionCorrectness = [](
    AdjacencyList&& adjacencyList, 
    const std::shared_ptr<NodeType>& comparisonTreePtr
  ) {
    auto treePtr = makeTree(adjacencyList);
    bool test = (*treePtr == *comparisonTreePtr);
    BOOST_CHECK(test);
    if(!test) {
      std::cout << "Expected tree: " << std::endl;
      std::cout << comparisonTreePtr << std::endl;
      std::cout << "Got: " << std::endl;
      std::cout << treePtr << std::endl << std::endl;
    }
  };

  // triangle
  testExpansionCorrectness(
    AdjacencyList(
      makeRepeatedElementCollection(Delib::ElementType::H, 3),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 2}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {2ul}),
      Tree::nodePtr(2ul, {1ul})
    })
  );

  // square
  testExpansionCorrectness(
    AdjacencyList(
      makeRepeatedElementCollection(Delib::ElementType::H, 4),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 3}, BondType::Single},
        {{2, 3}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {3ul}),
      Tree::nodePtr(2ul, {3ul})
    })
  );

  // pentangle
  testExpansionCorrectness(
    AdjacencyList(
      makeRepeatedElementCollection(Delib::ElementType::H, 5),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 3}, BondType::Single},
        {{2, 4}, BondType::Single},
        {{3, 4}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {
        Tree::nodePtr(3ul, {4ul})
      }),
      Tree::nodePtr(2ul, {
        Tree::nodePtr(4ul, {3ul})
      })
    })
  );

  // spiro for good measure
  testExpansionCorrectness(
    AdjacencyList(
      makeRepeatedElementCollection(Delib::ElementType::H, 7),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 3}, BondType::Single},
        {{2, 3}, BondType::Single},
        {{3, 4}, BondType::Single},
        {{3, 5}, BondType::Single},
        {{4, 6}, BondType::Single},
        {{5, 6}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {
        Tree::nodePtr(3ul, {
          Tree::nodePtr(4ul, {6ul}),
          Tree::nodePtr(5ul, {6ul})
        })
      }),
      Tree::nodePtr(2ul, {
        Tree::nodePtr(3ul, {
          Tree::nodePtr(4ul, {6ul}),
          Tree::nodePtr(5ul, {6ul})
        })
      })
    })
  );

  // fused triangles
  testExpansionCorrectness(
    AdjacencyList(
      makeRepeatedElementCollection(Delib::ElementType::H, 4),
      Edges({
        {{0, 1}, BondType::Single},
        {{0, 2}, BondType::Single},
        {{1, 2}, BondType::Single},
        {{1, 3}, BondType::Single},
        {{2, 3}, BondType::Single}
      })
    ),
    Tree::nodePtr(0ul, {
      Tree::nodePtr(1ul, {
        Tree::nodePtr(3ul),
        Tree::nodePtr(2ul, {3ul})
      }),
      Tree::nodePtr(2ul, {
        Tree::nodePtr(3ul),
        Tree::nodePtr(1ul, {3ul})
      }),
    })
  );

  { // order independence of tree generation on sequence in AdjacencyList
    Edges edges({
      {{0, 1}, BondType::Single},
      {{1, 2}, BondType::Single},
      {{1, 5}, BondType::Single},
      {{2, 3}, BondType::Single},
      {{3, 4}, BondType::Single},
      {{4, 5}, BondType::Single},
      {{4, 6}, BondType::Single},
      {{6, 7}, BondType::Single},
      {{6, 8}, BondType::Single}
    });

    auto prototypeTree = makeTree(
      AdjacencyList(
        makeRepeatedElementCollection(Delib::ElementType::H, 9),
        edges
      ), 
      0
    );

    std::vector<unsigned> sequence (edges.size());
    std::iota(
      sequence.begin(),
      sequence.end(),
      0
    );

    bool pass = true;
    do {
      // create new AdjacencyList and fill it with edges in the specified manner
      AdjacencyList adjacencies;

      for(unsigned i = 0; i < 9; i++) {
        adjacencies.addAtom(Delib::ElementType::H);
      }

      for(const auto& index : sequence) {
        auto edgesConstIterator = edges.begin();
        std::advance(edgesConstIterator, index);
        adjacencies.addBond(
          edgesConstIterator -> first.first,
          edgesConstIterator -> first.second,
          BondType::Single
        );
      }

      auto currentTree = makeTree(adjacencies, 0);

      if(*prototypeTree != *currentTree) {
        pass = false;
        break;
      }

    } while(std::next_permutation(sequence.begin(), sequence.end()));

    BOOST_CHECK(pass);
  }
}
