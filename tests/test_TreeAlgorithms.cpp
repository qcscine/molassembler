#include "BoostTestingHeader.h"
#include "TreeAlgorithms.h"
#include "AdjacencyListAlgorithms.h"

BOOST_AUTO_TEST_CASE( makeTreeTest ) {
  using namespace MoleculeManip;
  using namespace BasicTree;

  using NodeType = Node<AtomIndexType>;

  { 
    std::shared_ptr<NodeType> rootPtr = std::make_shared<NodeType>(
      0
    );

    //auto child = addChild<AtomIndexType>(rootPtr, 4);
    auto childPtr = rootPtr -> addChild(4);
    childPtr -> addChild(9);
    childPtr -> addChild(5);

    auto child2Ptr = rootPtr -> addChild(3);
    child2Ptr -> addChild(6);
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

    auto treePtr = AdjacencyListAlgorithms::makeTree(test);

    std::vector<AtomIndexType> BFSVisitSequence;
    auto BFSVisitor = [&BFSVisitSequence](const auto& nodePtr) {
      BFSVisitSequence.push_back(nodePtr -> key);
      return true;
    };
    BFSVisit(
      treePtr,
      BFSVisitor
    );

    BOOST_CHECK(BFSVisitSequence == std::vector<AtomIndexType>({
      0, 1, 2, 4, 3, 3, 5, 6, 7
    }));

    std::vector<AtomIndexType> DFSVisitSequence;
    auto DFSVisitor = [&DFSVisitSequence](const auto& nodePtr) {
      DFSVisitSequence.push_back(nodePtr -> key);
      return true;
    };

    DFSVisit(
      treePtr,
      DFSVisitor
    );

    BOOST_CHECK(DFSVisitSequence == std::vector<AtomIndexType>({
      0, 1, 4, 5, 7, 6, 3, 2, 3
    }));
  }
}
