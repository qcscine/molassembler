#define BOOST_TEST_MODULE RankingTreeTestsModule
#include "boost/test/unit_test.hpp"

#include "boost/graph/isomorphism.hpp"
#include "temple/Containers.h"
#include "temple/Stringify.h"

#include "molassembler/detail/StdlibTypeAlgorithms.h"
#include "molassembler/Cycles.h"
#include "molassembler/GraphAlgorithms.h"
#include "molassembler/IO.h"
#include "molassembler/RankingTree.h"
#include "molassembler/StereocenterList.h"

#include <random>
#include <fstream>

using namespace molassembler;

const std::string directoryPrefix = "test_files/ranking_tree_molecules/";

bool isStereocenter(
  const Molecule& molecule,
  const GraphType::edge_descriptor& e,
  const unsigned numPermutations,
  const boost::optional<unsigned>& assignment
) {
  auto stereocenterOption = molecule.getStereocenterList().option(e);

  if(!stereocenterOption) {
    std::cout << "No stereocenter on vertices " << temple::stringify(
      molecule.vertices(e)
    ) << "\n";
    return false;
  }

  if(stereocenterOption->numStereopermutations() != numPermutations) {
    std::cout << "Bond stereocenter on "
      << temple::stringify(
        molecule.vertices(e)
      )
      << " has " << stereocenterOption->numStereopermutations()
      << " stereopermutations, not " << numPermutations << "\n";
    return false;
  }

  if(assignment) {
    if(stereocenterOption->assigned() != assignment.value()) {
      std::cout << "Bond stereocenter on "
        << temple::stringify(
          molecule.vertices(e)
        )
        << " is assigned "
        << (
          stereocenterOption->assigned()
          ? std::to_string(stereocenterOption->assigned().value())
          : "u"
        )
        << ", not " << assignment.value() << "\n";
      return false;
    }
  }

  return true;
}

bool isStereocenter(
  const Molecule& molecule,
  AtomIndexType i,
  const unsigned numPermutations,
  const boost::optional<unsigned>& assignment
) {
  auto stereocenterOption = molecule.getStereocenterList().option(i);

  if(!stereocenterOption) {
    std::cout << "No stereocenter on atom index " << i << "\n";
    return false;
  }

  if(stereocenterOption->numStereopermutations() != numPermutations) {
    std::cout << "Atom stereocenter on " << i << " has "
      << stereocenterOption->numStereopermutations() << " stereopermutations, not "
      << numPermutations << "\n";
    return false;
  }

  if(assignment) {
    if(stereocenterOption->assigned() != assignment.value()) {
      std::cout << "Atom stereocenter on " << i << " is assigned "
        << (
          stereocenterOption->assigned()
          ? std::to_string(stereocenterOption->assigned().value())
          : "u"
        ) << ", not " << assignment.value() << "\n";
      return false;
    }
  }

  return true;
}

bool isStereogenic(
  const Molecule& molecule,
  AtomIndexType i
) {
  auto stereocenterOption = molecule.getStereocenterList().option(i);

  if(!stereocenterOption) {
    return false;
  }

  if(stereocenterOption->numStereopermutations() <= 1) {
    return false;
  }

  return true;
}

bool checkIsomorphicExpansion(
  const std::string& fileName,
  const AtomIndexType& expandOnIndex,
  const RankingTree::TreeGraphType& comparisonGraph
) {
  auto molecule = IO::read(
    directoryPrefix + fileName
  );

  auto expandedTree = RankingTree(
    molecule.getGraph(),
    molecule.getCycleData(),
    molecule.getStereocenterList(),
    molecule.dumpGraphviz(),
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
  auto molecule = IO::read(
    directoryPrefix + fileName
  );

  auto expandedTree = RankingTree(
    molecule.getGraph(),
    molecule.getCycleData(),
    molecule.getStereocenterList(),
    molecule.dumpGraphviz(),
    expandOnIndex,
    {},
    RankingTree::ExpansionOption::Full
  );

  std::ofstream dotFile(fileName + ".dot");
  dotFile << expandedTree.dumpGraphviz();
  dotFile.close();
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

BOOST_AUTO_TEST_CASE(TreeExpansionAndSequenceRuleOneTests) {
  using namespace std::string_literals;


  // Basic tests

  /* P-92.2.1.1.2 Spheres I and II */
  auto exampleOne = IO::read(
    directoryPrefix + "2R-2-chloropropan-1-ol.mol"s
  );

  auto exampleOneExpanded = RankingTree(
    exampleOne.getGraph(),
    exampleOne.getCycleData(),
    exampleOne.getStereocenterList(),
    exampleOne.dumpGraphviz(),
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


  auto exampleTwo = IO::read(
    directoryPrefix + "2S-23-dichloropropan-1-ol.mol"
  );

  auto exampleTwoExpanded = RankingTree(
    exampleTwo.getGraph(),
    exampleTwo.getCycleData(),
    exampleTwo.getStereocenterList(),
    exampleTwo.dumpGraphviz(),
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
  auto exampleThree = IO::read(
    directoryPrefix + "1S5R-bicyclo-3-1-0-hex-2-ene.mol"
  );

  auto exampleThreeExpanded = RankingTree(
    exampleThree.getGraph(),
    exampleThree.getCycleData(),
    exampleThree.getStereocenterList(),
    exampleThree.dumpGraphviz(),
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
    << temple::condenseIterable(
      temple::map(
        exampleThreeRanked,
        [](const auto& set) -> std::string {
          return temple::condenseIterable(set);
        }
      )
    )
  );

  auto exampleThreeExpandedAgain = RankingTree(
    exampleThree.getGraph(),
    exampleThree.getCycleData(),
    exampleThree.getStereocenterList(),
    exampleThree.dumpGraphviz(),
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
std::string condenseSets(const std::vector<std::vector<T>>& sets) {
  return temple::condenseIterable(
    temple::map(
      sets,
      [](const auto& set) -> std::string {
        return "{"s + temple::condenseIterable(set) + "}"s;
      }
    )
  );
}

BOOST_AUTO_TEST_CASE(sequenceRuleThreeTests) {
  // P-92.4.2.1 Example 1 (Z before E)
  auto ZEDifference = IO::read(
    directoryPrefix + "2Z5S7E-nona-2,7-dien-5-ol.mol"s
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(ZEDifference, 0, 2, 0),
    "Stereocenter at C0 in 2Z5S7E-nona-2,7-dien-5-ol is not S"
  );

  // P-92.4.2.2 Example 1 (Z before E in aux. stereocenters, splitting)
  auto EECyclobutane = IO::read(
    directoryPrefix + "1E3E-1,3-difluoromethylidenecyclobutane.mol"
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(EECyclobutane, EECyclobutane.edge(0, 3), 2, 0)
    && isStereocenter(EECyclobutane, EECyclobutane.edge(5, 6), 2, 0),
    "1E3E-1,3-difluoromethylidenecyclobutane double bonds aren't E"
  );

  // P-92.4.2.2 Example 2 (stereogenic before non-stereogenic)
  auto inTreeNstgDB = IO::read(
    directoryPrefix
    + "(2Z5Z7R8Z11Z)-9-(2Z-but-2-en-1-yl)-5-(2E-but-2-en-1-yl)trideca-2,5,8,11-tetraen-7-ol.mol"s
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(inTreeNstgDB, 0, 2, 1u),
    "(2Z5Z7R8Z11Z)-9-(2Z-but-2-en-1-yl)-5-(2E-but-2-en-1-yl)trideca-2,5,8,11-tetraen-7-ol "
    "difference between non-stereogenic auxiliary stereocenter and assigned "
    "stereocenter isn't recognized!"
  );
}

BOOST_AUTO_TEST_CASE(sequenceRuleFourTests) {
  /* TODO
   * is it necessary to add a test to ensure full partial ordering?
   * stereogenic > pseudostereogenic > non-stereogenic
   *
   * currently no differentiation between stereogenic and pseudostereogenic
   */

  // (4A) P-92.5.1 Example (stereogenic before non-stereogenic)
  auto pseudoOverNonstg = IO::read(
    directoryPrefix
    + "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol"s
  );

  BOOST_CHECK_MESSAGE(
    !isStereogenic(pseudoOverNonstg, 10),
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "branch with R-R aux. stereocenters not non-stereogenic"
  );

  BOOST_CHECK_MESSAGE(
    isStereogenic(pseudoOverNonstg, 1),
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "branch with R-S aux. stereocenters not stereogenic"
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(pseudoOverNonstg, 0, 2, 0),
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "sequence rule 4A does not recognize stereogenic over non-stereogenic, 3 as S"
  );

  // (4B) P-92.5.2.2 Example 1 (single chain pairing, ordering and reference selection)
  auto simpleLikeUnlike = IO::read(
    directoryPrefix + "(2R,3R,4R,5S,6R)-2,3,4,5,6-pentachloroheptanedioic-acid.mol"s
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(simpleLikeUnlike, 10, 2, 1u),
    "(2R,3R,4R,5S,6R)-2,3,4,5,6-pentachloroheptanedioic-acid central carbon does "
    " not register as a stereocenter and/or isn't assigned as R"
  );

  // (4B) P-92.5.2.2 Example 3 (single-chain pairing, cycle splitting)
  auto lAlphaLindane = IO::read(
    directoryPrefix + "l-alpha-lindane.mol"s
  );

  BOOST_CHECK_MESSAGE(
    (
      temple::all_of(
        std::vector<AtomIndexType> {6, 7, 8, 9, 10, 11},
        [&](const auto carbonIndex) -> bool {
          return isStereogenic(lAlphaLindane, carbonIndex);
        }
      )
    ),
    "Not all L-alpha-lindane carbon atoms not recognized as stereocenters!"
  );

  // (4B) P-92.5.2.2 Example 4 (multiple-chain stereocenter ranking)
  auto oxyNitroDiffBranches = IO::read(
    directoryPrefix + "(2R,3S,6R,9R,10S)-6-chloro-5-(1R,2S)-1,2-dihydroxypropoxy-7-(1S,2S)-1,2-dihydroxypropoxy-4,8-dioxa-5,7-diazaundecande-2,3,9,10-tetrol.mol"s
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(oxyNitroDiffBranches, 0, 2, 1u),
    "(2R,3S,6R,9R,10S)-6-chloro-5-(1R,2S)-1,2-dihydroxypropoxy-7-(1S,2S)-1,2-dihydroxypropoxy-4,8-dioxa-5,7-diazaundecande-2,3,9,10-tetrol central carbon not recognized as R"
  );

  // (4B) P-92.5.2.2 Example 5 (multiple-chain stereocenter ranking)
  auto groupingDifferences = IO::read(
    directoryPrefix + "(2R,3R,5R,7R,8R)-4.4-bis(2S,3R-3-chlorobutan-2-yl)-6,6-bis(2S,4S-3-chlorobutan-2-yl)-2,8-dichloro-3,7-dimethylnonan-5-ol.mol"s
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(groupingDifferences, 0, 2, 1u),
    "The central carbon in (2R,3R,5R,7R,8R)-4.4-bis(2S,3R-3-chlorobutan-2-yl)-6,6-bis(2S,4S-3-chlorobutan-2-yl)-2,8-dichloro-3,7-dimethylnonan-5-ol is not recognized as R"
  );

  // (4B) P-92.5.2.2 Example 6 (number of reference descriptors)
  auto numReferenceDescriptors = IO::read(
    directoryPrefix + "2R-2-bis(1R)-1-hydroxyethylamino-2-(1R)-1-hydroxyethyl(1S)-1-hydroxyethylaminoacetic-acid.mol"
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(numReferenceDescriptors, 0, 2, 1u),
    "The central carbon in 2R-2-bis(1R)-1-hydroxyethylamino-2-(1R)-1-hydroxyethyl(1S)-1-hydroxyethylaminoacetic-acid is not recognized as R"
  );
}

BOOST_AUTO_TEST_CASE(sequenceRuleFiveTests) {
  // (4C) P-92.5.3 Example r/s leads to R difference
  auto rsDifference = IO::read(
    directoryPrefix + "(2R,3r,4R,5s,6R)-2,6-dichloro-3,5-bis(1S-1-chloroethyl)heptan-4-ol.mol"
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(rsDifference, 0, 2, 1u),
    "The central carbon in (2R,3r,4R,5s,6R)-2,6-dichloro-3,5-bis(1S-1-chloroethyl)heptan-4-ol is not recognized as R"
  );

  // (5) P-92.6 Example 1 simple R/S difference leads to r
  auto pseudo = IO::read(
    directoryPrefix + "(2R,3r,4S)-pentane-2,3,4-trithiol.mol"
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(pseudo, 0, 2, 1u),
    "The central carbon in (2R,3r,4S)-pentane-2,3,4-trithiol is not recognized as R"
  );

  // (5) P-92.6 Example 2 cyclobutane splitting
  auto cyclobutane = IO::read(
    directoryPrefix + "(1r,3r)-cyclobutane-1,3-diol.mol"
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(cyclobutane, 2, 2, 1u)
    && isStereocenter(cyclobutane, 3, 2, 1u),
    "The chiral carbons in (1r,3r)-cyclobutane-1,3-diol aren't properly recognized"
  );

  // (5) P-92.6 Example 5 double bond ranking
  auto pseudoDB = IO::read(
    directoryPrefix + "(2E,4R)-4-chloro-3-(1S-1-chloroethyl)pent-2-ene.mol"
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(pseudoDB, 0, 2, 0u),
    "Double bond in (2E,4R)-4-chloro-3-(1S-1-chloroethyl)pent-2-ene isn't E"
  );

  // (5) P-92.6 Example 6
  auto fourDoesNothing = IO::read(
    directoryPrefix + "1s-1-(1R,2R-1,2-dichloropropyl-1S,2R-1,2-dichloropropylamino)1-(1R,2S-1,2-dichloropropyl-1S,2S-1,2-dichloropropylamino)methan-1-ol.mol"
  );

  BOOST_CHECK_MESSAGE(
    isStereocenter(fourDoesNothing, 0, 2, 0u),
    "The central stereocenter in 1s-1-(1R,2R-1,2-dichloropropyl-1S,2R-1,2-dichloropropylamino)1-(1R,2S-1,2-dichloropropyl-1S,2S-1,2-dichloropropylamino)methan-1-ol isn't recognized as S"
  );
}
