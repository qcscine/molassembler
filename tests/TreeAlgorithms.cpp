#define BOOST_TEST_MODULE TreeAlgorithmsTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "TreeAlgorithms.h"
#include "AdjacencyListAlgorithms.h"
#include "RepeatedElementCollection.h"

BOOST_AUTO_TEST_CASE( makeTreeTest ) {
  using namespace MoleculeManip;
  using namespace Tree;

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

    auto treePtr = AdjacencyListAlgorithms::makeTree(test);

    std::vector<AtomIndexType> BFSVisitSequence;
    auto BFSVisitor = [&BFSVisitSequence](const auto& nodePtr) {
      BFSVisitSequence.push_back(nodePtr -> key);
      return true;
    };

    TreeAlgorithms::BFSVisit(
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

    TreeAlgorithms::DFSVisit(
      treePtr,
      DFSVisitor
    );

    BOOST_CHECK(DFSVisitSequence == std::vector<AtomIndexType>({
      0, 1, 4, 5, 7, 6, 3, 2, 3
    }));
  }
  { // BFS, DFS with depth limit testing
    AdjacencyList test(
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

    auto treePtr = AdjacencyListAlgorithms::makeTree(test);

    std::vector<AtomIndexType> BFSVisitSequence;
    auto BFSVisitor = [&BFSVisitSequence](const auto& nodePtr) {
      BFSVisitSequence.push_back(nodePtr -> key);
      return true;
    };

    TreeAlgorithms::BFSVisit(
      treePtr,
      BFSVisitor,
      3
    );

    BOOST_CHECK(BFSVisitSequence == std::vector<AtomIndexType>({
      0, 1, 2, 4, 3, 3, 5
    }));

    std::vector<AtomIndexType> DFSVisitSequence;
    auto DFSVisitor = [&DFSVisitSequence](const auto& nodePtr) {
      DFSVisitSequence.push_back(nodePtr -> key);
      return true;
    };

    TreeAlgorithms::DFSVisit(
      treePtr,
      DFSVisitor,
      3
    );

    BOOST_CHECK(DFSVisitSequence == std::vector<AtomIndexType>({
      0, 1, 4, 5, 3, 2, 3
    }));
  }
  {
    AdjacencyList test(
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

    auto treePtr = AdjacencyListAlgorithms::makeTree(test);

    struct BFSVisitor {
      std::vector<AtomIndexType> visitSequence;

      bool operator() (
        const std::shared_ptr<NodeType>& nodePtr,
        const unsigned& depth __attribute__ ((unused))
      ) {
        visitSequence.push_back(
          nodePtr -> key
        );
        return true;
      }
    };

    BFSVisitor visitor;

    TreeAlgorithms::BFSVisit(
      treePtr,
      visitor,
      3
    );

    BOOST_CHECK(visitor.visitSequence == std::vector<AtomIndexType>({
      0, 1, 2, 4, 3, 3, 5
    }));
  }
}
