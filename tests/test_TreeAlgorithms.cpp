#include "BoostTestingHeader.h"
#include "TreeAlgorithms.h"

BOOST_AUTO_TEST_CASE( makeTreeTest ) {
  using namespace MoleculeManip;
  using namespace BasicTree;

  { // output
    std::shared_ptr<NodeType> rootPtr = std::make_shared<NodeType>(
      0
    );

    //auto child = addChild<AtomIndexType>(rootPtr, 4);
    auto childPtr = rootPtr -> addChild(4);
    childPtr -> addChild(9);
    childPtr -> addChild(5);

    auto child2Ptr = rootPtr -> addChild(3);
    child2Ptr -> addChild(6);

    std::string expectedString = "digraph tree {\n  \"0\" -> \"4\";\n"
      "  \"0\" -> \"3\";\n  \"4\" -> \"9\";\n  \"4\" -> \"5\";\n"
      "  \"3\" -> \"6\";\n}\n";

    BOOST_CHECK(rootPtr -> toString() == expectedString);
  }
  { // expansion
    AdjacencyList test(
      EdgeList({
        Edge(0, 1, BondType::Single),
        Edge(0, 2, BondType::Single),
        Edge(1, 2, BondType::Single)
      })
    );

    std::string expectedString = "digraph tree {\n  \"0\" -> \"1\";\n"
      "  \"1\" -> \"2\";\n  \"2\" -> \"0\";\n}\n";

    auto treeStruct = makeTree(test);

    BOOST_CHECK(treeStruct.rootPtr -> toString() == expectedString);
  }
}
