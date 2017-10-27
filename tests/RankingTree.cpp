#define BOOST_TEST_MODULE RankingTreeTestsModule
#define BOOST_TEST_DYN_LINK
#include "boost/test/unit_test.hpp"

#include "boost/graph/isomorphism.hpp"
#include "template_magic/Containers.h"

#include "GraphAlgorithms.h"
#include "IO.h"
#include "RankingTree.h"
#include "RepeatedElementCollection.h"
#include "StdlibTypeAlgorithms.h"

#include <random>
#include <fstream>

/* TODO
 * - add tests from every sequence rule that is applied and see if the resulting
 *   ranking is correct!
 */

using namespace MoleculeManip;

bool checkIsomorphicExpansion(
  const std::string& fileName,
  const AtomIndexType& expandOnIndex,
  const RankingTree::TreeGraphType& comparisonGraph
) {
  IO::MOLFileHandler molHandler;
  auto molecule = molHandler.readSingle(
    "../tests/mol_files/ranking_tree_molecules/"s
    + fileName
  );

  auto expandedTree = RankingTree(molecule, expandOnIndex);

  return boost::graph::isomorphism(
    expandedTree.getGraph(),
    comparisonGraph
  );
}

void writeExpandedTree(
  const std::string& fileName,
  const AtomIndexType& expandOnIndex
) {
  IO::MOLFileHandler molHandler;
  auto molecule = molHandler.readSingle(
    "../tests/mol_files/ranking_tree_molecules/"s
    + fileName
  );

  auto expandedTree = RankingTree(molecule, expandOnIndex);

  std::ofstream dotFile(fileName + ".dot");
  dotFile << expandedTree.dumpGraphviz();
  dotFile.close();
}

void runBFSExplainer(
  const std::string& fileName,
  const AtomIndexType& startOn
) {
  IO::MOLFileHandler molHandler;
  auto molecule = molHandler.readSingle(
    "../tests/mol_files/ranking_tree_molecules/"s
    + fileName
  );

  GraphAlgorithms::BFSVisitors::ShowAllEventsBFSVisitor visitor;

  using ColorMapBase = std::map<
    AtomIndexType,
    boost::default_color_type
  >;

  ColorMapBase colorMap;
  boost::associative_property_map<ColorMapBase> propColorMap(colorMap);
  boost::queue<GraphType::vertex_descriptor> Q;

  boost::breadth_first_visit(
    // The graph to operate on
    molecule.getGraph(),
    // The vertex to start with
    startOn,
    // A queue object to store vertex_descriptors
    Q,
    // The visitor to use
    visitor,
    // A map to store color (state)
    propColorMap
  );
}

RankingTree::TreeGraphType makeTree(
  const std::vector<
    std::pair<unsigned, unsigned>
  >& edges
) {
  RankingTree::TreeGraphType tree;

  auto addNodeIfNotExists = [&](const unsigned& nodeIndex) {
    while(nodeIndex >= boost::num_vertices(tree)) {
      boost::add_vertex(tree);
    }
  };

  for(const auto& pair : edges) {
    addNodeIfNotExists(pair.first);
    addNodeIfNotExists(pair.second);

    boost::add_edge(pair.first, pair.second, tree);
  }

  return tree;
}

BOOST_AUTO_TEST_CASE(trivialExamples) {
  using namespace std::string_literals;

  IO::MOLFileHandler molHandler;
  // TODO Continue here

}

BOOST_AUTO_TEST_CASE(IUPAC2013Examples) {
  using namespace std::string_literals;

  IO::MOLFileHandler molHandler;
  std::string directoryPrefix = "../tests/mol_files/ranking_tree_molecules/"s;

  // Basic tests

  /* P-92.2.1.1.2 Spheres I and II */
  auto exampleOne = molHandler.readSingle(
    directoryPrefix
    + "2R-2-chloropropan-1-ol.mol"s
  );

  auto exampleOneExpanded = RankingTree(exampleOne, 2);

  BOOST_CHECK(
    boost::graph::isomorphism(
      exampleOneExpanded.getGraph(),
      makeTree({
        {0, 1},
        {0, 2},
        {0, 3},
        {0, 4},
        {3, 5},
        {3, 6},
        {3, 7},
        {4, 8},
        {4, 9},
        {4, 10},
        {10, 11}
      })
    )
  );
      

  auto exampleTwo = molHandler.readSingle(
    directoryPrefix
    + "2S-23-dichloropropan-1-ol.mol"
  );

  auto exampleTwoExpanded = RankingTree(exampleTwo, 3);

  BOOST_CHECK(
    boost::graph::isomorphism(
      exampleTwoExpanded.getGraph(),
      makeTree({
        {0, 1}, // Cl
        {0, 2}, // CH2OH
        {0, 3}, // CH2Cl
        {0, 4},
        {2, 5}, // CO
        {2, 6},
        {2, 7},
        {5, 8}, // OH
        {3, 9}, // Cl
        {3, 10},
        {3, 11},
        {4, 12},
        {4, 13},
        {4, 14}
      })
    )
  );

  // P. 92.2.2 Sequence subrule 1b: Priority due to duplicate atoms
  // Cycle and multiple-bond splitting
  auto exampleThree = molHandler.readSingle(
    directoryPrefix
    + "1S5R-bicyclo-3-1-0-hex-2-ene.mol"
  );

  auto exampleThreeExpanded = RankingTree(exampleThree, 0);

  BOOST_CHECK(
    boost::graph::isomorphism(
      exampleThreeExpanded.getGraph(),
      makeTree({
        { 0,   1}, // 0 -> 1
        { 1,   2}, //      └> 4
        { 2,   3}, //      |  └> (5)
        { 2,   4}, //      |  └> 12
        { 2,   5}, //      |  └> 5
        { 5,   6}, //      |     └> (4)
        { 5,   7}, //      |     └> 13
        { 5,   8}, //      |     └> 3
        { 8,   9}, //      |        └> 10
        { 8,  10}, //      |        └> 11
        { 8,  11}, //      |        └> (0)
        { 1,  12}, //      └> 7
        { 1,  13}, //      └> 2
        {13,  14}, //         └> (0)
        {13,  15}, //         └> 8
        {13,  16}, //         └> 9
        { 0,  17}, // 0 -> 2
        {17,  18}, //      └> 8
        {17,  19}, //      └> 9
        {17,  20}, //      └> 1
        {20,  21}, //         └> (0)
        {20,  22}, //         └> 7
        {20,  23}, //         └> 4
        {23,  24}, //            └> (5)
        {23,  25}, //            └> 12
        {23,  26}, //            └> 5
        {26,  27}, //               └> (4)
        {26,  28}, //               └> 13
        {26,  29}, //               └> 3
        {29,  30}, //                  └> 10
        {29,  31}, //                  └> 11
        {29,  32}, //                  └> (0)
        { 0,  33}, // 0 -> 3
        {33,  34}, //      └> 10
        {33,  35}, //      └> 11
        {33,  36}, //      └> 5
        {36,  37}, //         └> (4)
        {36,  38}, //         └> 13
        {36,  39}, //         └> 4
        {39,  40}, //            └> (5)
        {39,  41}, //            └> 12
        {39,  42}, //            └> 1
        {42,  43}, //               └> (0)
        {42,  44}, //               └> 7
        {42,  45}, //               └> 2
        {45,  46}, //                  └> (0)
        {45,  47}, //                  └> 8
        {45,  48}, //                  └> 9
        { 0,  49}, // 0 -> 6
      })
    )
  );

  auto exampleThreeRanked = exampleThreeExpanded.rank();

  BOOST_CHECK_MESSAGE(
    (exampleThreeRanked == std::vector<
      std::vector<unsigned long>
    > { {6}, {3}, {2}, {1} }),
    "Example three expanded is not {{6}, {3}, {2}, {1}}, but: "
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        exampleThreeRanked,
        [](const auto& set) -> std::string {
          return TemplateMagic::condenseIterable(set);
        }
      )
    )
  );

  auto exampleThreeExpandedAgain = RankingTree(exampleThree, 1);

  BOOST_CHECK((
    exampleThreeExpandedAgain.rank() == std::vector<
      std::vector<unsigned long>
    > { {7}, {4}, {2}, {0} }
  ));

  // Write out the ranked result
  /*std::cout << TemplateMagic::condenseIterable(
    TemplateMagic::map(
      reRanked,
      [](const auto& set) -> std::string {
        return "{"s + TemplateMagic::condenseIterable(set) + "}"s;
      }
    )
  ) << std::endl;*/
}
