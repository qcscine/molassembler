#include "BoostTestingHeader.h"

#include "AdjacencyListAlgorithms.h"
#include "StdlibTypeAlgorithms.h"


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
      EdgeList({
        Edge(0, 1, BondType::Single),
        Edge(0, 2, BondType::Single),
        Edge(1, 2, BondType::Single)
      })
    );

    std::string expectedString = "digraph tree {\n"
      "  aa[label=\"0\"];\n"
      "  ab[label=\"1\"];\n"
      "  ac[label=\"2\"];\n"
      "  ad[label=\"2\"];\n"
      "  aa -> ab;\n"
      "  aa -> ac;\n"
      "  ab -> ad;\n"
      "}\n";

    auto treePtr = makeTree(test);

    BOOST_CHECK(treePtr -> toString() == expectedString);
  }
  { // BFS, DFS testing
    AdjacencyList test(
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

    auto treePtr = makeTree(test, 0);

    std::string expectedString = "digraph tree {\n"
      "  aa[label=\"0\"];\n"
      "  ab[label=\"1\"];\n"
      "  ac[label=\"2\"];\n"
      "  ad[label=\"4\"];\n"
      "  ae[label=\"3\"];\n"
      "  af[label=\"3\"];\n"
      "  ag[label=\"5\"];\n"
      "  ah[label=\"6\"];\n"
      "  ai[label=\"7\"];\n"
      "  aa -> ab;\n"
      "  ab -> ac;\n"
      "  ab -> ad;\n"
      "  ac -> ae;\n"
      "  ad -> af;\n"
      "  ad -> ag;\n"
      "  ag -> ah;\n"
      "  ag -> ai;\n"
      "}\n";

    BOOST_CHECK(treePtr -> toString() == expectedString);
  }
}
