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

  auto expandedTree = RankingTree(
    molecule,
    expandOnIndex,
    {},
    RankingTree::ExpansionOption::Full
  );

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

  auto expandedTree = RankingTree(
    molecule,
    expandOnIndex,
    {},
    RankingTree::ExpansionOption::Full
  );

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

const std::string directoryPrefix = "../tests/mol_files/ranking_tree_molecules/"s;

BOOST_AUTO_TEST_CASE(TreeExpansionAndSequenceRuleOneTests) {
  using namespace std::string_literals;

  IO::MOLFileHandler molHandler;

  // Basic tests

  /* P-92.2.1.1.2 Spheres I and II */
  auto exampleOne = molHandler.readSingle(
    directoryPrefix + "2R-2-chloropropan-1-ol.mol"s
  );

  auto exampleOneExpanded = RankingTree(
    exampleOne,
    2,
    {},
    RankingTree::ExpansionOption::Full
  );

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
    directoryPrefix + "2S-23-dichloropropan-1-ol.mol"
  );

  auto exampleTwoExpanded = RankingTree(
    exampleTwo,
    3,
    {},
    RankingTree::ExpansionOption::Full
  );

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
    directoryPrefix + "1S5R-bicyclo-3-1-0-hex-2-ene.mol"
  );

  auto exampleThreeExpanded = RankingTree(
    exampleThree,
    0,
    {},
    RankingTree::ExpansionOption::Full
  );

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

  auto exampleThreeRanked = exampleThreeExpanded.getRanked();

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

  auto exampleThreeExpandedAgain = RankingTree(
    exampleThree,
    1,
    {},
    RankingTree::ExpansionOption::Full
  );

  BOOST_CHECK((
    exampleThreeExpandedAgain.getRanked() == std::vector<
      std::vector<unsigned long>
    > { {7}, {4}, {2}, {0} }
  ));
}

template<typename T>
std::string condenseSets (const std::vector<std::vector<T>>& sets) {
  return TemplateMagic::condenseIterable(
    TemplateMagic::map(
      sets,
      [](const auto& set) -> std::string {
        return "{"s + TemplateMagic::condenseIterable(set) + "}"s;
      }
    )
  );
}

BOOST_AUTO_TEST_CASE(sequenceRuleThreeTests) {
  IO::MOLFileHandler molHandler;

  // P-92.4.2.1 Example 1 (Z before E)
  auto ZEDifference = molHandler.readSingle(
    directoryPrefix + "2Z5S7E-nona-2,7-dien-5-ol.mol"s
  );

  const auto& stereocenters = ZEDifference.getStereocenterList();

  BOOST_CHECK_MESSAGE(
    stereocenters.involving(0)
    && stereocenters.at(0)->numAssignments() == 2,
    "Stereocenter at C0 in 2Z5S7E-nona-2,7-dien-5-ol is not found "
    "or is not determined as having two assignments."
  );

  BOOST_CHECK_MESSAGE(
    stereocenters.at(0)->assigned() == 0u,
    "Stereocenter at C0 in 2Z5S7E-nona-2,7-dien-5-ol is not S"
  );

  /* Re-expand (assignment of auxiliary EZStereocenters has to occur from
   * molecule information
   */
  auto reExpanded = RankingTree(ZEDifference, 0);

  BOOST_CHECK_MESSAGE(
    reExpanded.getRanked().size() == 4,
    "Re-expanding 2Z5S7E-nona-2,7-dien-5-ol at 0 does not yield a difference "
    " between the Z and E branches. Perhaps the auxiliary stereocenters aren't "
    " properly instantiated from molecular graph information in sequence rule 3 "
    " prep? The ranked sets are " << condenseSets(reExpanded.getRanked())
  );

  // P-92.4.2.2 Example 1 (Z before E in aux. stereocenters, splitting)
  auto EECyclobutane = molHandler.readSingle(
    directoryPrefix + "1E3E-1,3-difluoromethylidenecyclobutane.mol"
  );

  BOOST_CHECK_MESSAGE(
    EECyclobutane.getStereocenterList().involving(0)
    && EECyclobutane.getStereocenterList().at(0)->numAssignments() == 2
    && EECyclobutane.getStereocenterList().involving(5)
    && EECyclobutane.getStereocenterList().at(5)->numAssignments() == 2,
    "1E3E-1,3-difluoromethylidenecyclobutane differences between branches "
    " when breaking the cyclobutane at the EZStereocenters don't register!"
  );

  BOOST_CHECK_MESSAGE(
    EECyclobutane.getStereocenterList().at(0)->assigned() == 0u
    && EECyclobutane.getStereocenterList().at(5)->assigned() == 0u,
    "1E3E-1,3-difluoromethylidenecyclobutane double bonds aren't E"
  );

  // P-92.4.2.2 Example 2 (stereogenic before non-stereogenic)
  auto inTreeNstgDB = molHandler.readSingle(
    directoryPrefix 
    + "(2Z5Z7R8Z11Z)-9-(2Z-but-2-en-1-yl)-5-(2E-but-2-en-1-yl)trideca-2,5,8,11-tetraen-7-ol.mol"s
  );

  BOOST_CHECK_MESSAGE(
    inTreeNstgDB.getStereocenterList().involving(0)
    && inTreeNstgDB.getStereocenterList().at(0)->numAssignments() == 2
    && inTreeNstgDB.getStereocenterList().at(0)->assigned() == 1u,
    "(2Z5Z7R8Z11Z)-9-(2Z-but-2-en-1-yl)-5-(2E-but-2-en-1-yl)trideca-2,5,8,11-tetraen-7-ol "
    "difference between non-stereogenic auxiliary stereocenter and assigned "
    "stereocenter isn't recognized!"
  );
}

BOOST_AUTO_TEST_CASE(sequenceRuleFourTests) {
  IO::MOLFileHandler molHandler;

  /* TODO
   * is it necessary to add a test to ensure full partial ordering?
   * stereogenic > pseudostereogenic > non-stereogenic
   *
   * currently no differentiation between stereogenic and pseudostereogenic
   */

  // (4A) P-92.5.1 Example (stereogenic before non-stereogenic)
  auto pseudoOverNonstg = molHandler.readSingle(
    directoryPrefix 
    + "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol"s
  );

  const auto& pseudoOverNonstgStereocenters = pseudoOverNonstg.getStereocenterList();

  BOOST_CHECK_MESSAGE(
    pseudoOverNonstgStereocenters.involving(10)
    && pseudoOverNonstgStereocenters.at(10)->numAssignments() == 1,
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "branch with R-R aux. stereocenters not non-stereogenic"
  );

  BOOST_CHECK_MESSAGE(
    pseudoOverNonstgStereocenters.involving(1)
    && pseudoOverNonstgStereocenters.at(1)->numAssignments() == 2,
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "branch with R-S aux. stereocenters not stereogenic"
  );

  BOOST_CHECK_MESSAGE(
    pseudoOverNonstgStereocenters.involving(0)
    && pseudoOverNonstgStereocenters.at(0)->numAssignments() == 2,
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "sequence rule 4A does not recognize stereogenic over non-stereogenic"
  );

  // (4B) P-92.5.2.2 Example 1 (single chain pairing, ordering and reference selection)
  auto simpleLikeUnlike = molHandler.readSingle(
    directoryPrefix + "(2R,3R,4R,5S,6R)-2,3,4,5,6-pentachloroheptanedioic-acid.mol"s
  );

  const auto& simpleLikeUnlikeStereocenters = simpleLikeUnlike.getStereocenterList();

  BOOST_CHECK_MESSAGE(
    simpleLikeUnlikeStereocenters.involving(10)
    && simpleLikeUnlikeStereocenters.at(10)->numAssignments() == 2
    && simpleLikeUnlikeStereocenters.at(10)->assigned() == 1u,
    "(2R,3R,4R,5S,6R)-2,3,4,5,6-pentachloroheptanedioic-acid central carbon does "
    " not register as a stereocenter and/or isn't assigned as R"
  );

  // (4B) P-92.5.2.2 Example 3 (single-chain pairing, cycle splitting)
  auto lAlphaLindane = molHandler.readSingle(
    directoryPrefix + "l-alpha-lindane.mol"s
  );

  const auto& lAlphaLindaneStereocenters = lAlphaLindane.getStereocenterList();

  BOOST_CHECK_MESSAGE(
    (
      TemplateMagic::all_of(
        TemplateMagic::map(
          std::vector<AtomIndexType> {6, 7, 8, 9, 10, 11},
          [&](const auto& carbonIndex) -> bool {
            return (
              lAlphaLindaneStereocenters.involving(carbonIndex)
              && lAlphaLindaneStereocenters.at(carbonIndex)->numAssignments() == 2
            );
          }
        )
      )
    ),
    "Not all L-alpha-lindane carbon atoms not recognized as stereocenters!"
  );

  // (4B) P-92.5.2.2 Example 4 (multiple-chain stereocenter ranking)
  auto oxyNitroDiffBranches = molHandler.readSingle(
    directoryPrefix + "(2R,3S,6R,9R,10S)-6-chloro-5-(1R,2S)-1,2-dihydroxypropoxy-7-(1S,2S)-1,2-dihydroxypropoxy-4,8-dioxa-5,7-diazaundecande-2,3,9,10-tetrol.mol"s
  );

  const auto& oxyNitroDiffBranchesStereocenters = oxyNitroDiffBranches.getStereocenterList();

  BOOST_CHECK_MESSAGE(
    oxyNitroDiffBranchesStereocenters.involving(0)
    && oxyNitroDiffBranchesStereocenters.at(0)->numAssignments() == 2
    && oxyNitroDiffBranchesStereocenters.at(0)->assigned() == 1u,
    "(2R,3S,6R,9R,10S)-6-chloro-5-(1R,2S)-1,2-dihydroxypropoxy-7-(1S,2S)-1,2-dihydroxypropoxy-4,8-dioxa-5,7-diazaundecande-2,3,9,10-tetrol central carbon not recognized as R"
  );

  // (4B) P-92.5.2.2 Example 5 (multiple-chain stereocenter ranking)
  auto groupingDifferences = molHandler.readSingle(
    directoryPrefix + "(2R,3R,5R,7R,8R)-4.4-bis(2S,3R-3-chlorobutan-2-yl)-6,6-bis(2S,4S-3-chlorobutan-2-yl)-2,8-dichloro-3,7-dimethylnonan-5-ol.mol"s
  );

  const auto& groupingDifferencesStereocenters = groupingDifferences.getStereocenterList();

  /* TODO
   * - R/S assignment of C0 is incorrect! When choosing representative
   *   stereocenters, R is chosen for both branches since sequence rules 4/5 are
   *   applied for both downwards from their respective junctions when comparing
   *   stereocenter groups. Would need to limit auxiliaryRank application depth
   *   in those situations to the depth prior to encountering the junction call
   *   subjects...
   */
  BOOST_CHECK_MESSAGE(
    groupingDifferencesStereocenters.involving(0)
    && groupingDifferencesStereocenters.at(0) -> numAssignments() == 2
    && groupingDifferencesStereocenters.at(0) -> assigned() == 1u,
    "The central carbon in (2R,3R,5R,7R,8R)-4.4-bis(2S,3R-3-chlorobutan-2-yl)-6,6-bis(2S,4S-3-chlorobutan-2-yl)-2,8-dichloro-3,7-dimethylnonan-5-ol is not recognized as R"
  );
}
