/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/test/unit_test.hpp"

#include "temple/constexpr/Bitmask.h"
#include "temple/Functional.h"
#include "temple/Stringify.h"

#include "molassembler/Cycles.h"
#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/IO.h"
#include "molassembler/Molecule/RankingTree.h"
#include "molassembler/StereopermutatorList.h"

#include <random>
#include <fstream>

using namespace Scine;
using namespace molassembler;

bool isBondStereopermutator(
  const Molecule& molecule,
  const BondIndex& e,
  const unsigned numPermutations,
  const boost::optional<unsigned>& assignment
) {
  auto stereopermutatorOption = molecule.stereopermutators().option(e);

  if(!stereopermutatorOption) {
    std::cout << "No stereopermutator on vertices " << e.first << " - " << e.second << "\n";

    for(const auto& stereopermutator : molecule.stereopermutators().bondStereopermutators()) {
      std::cout << "BondStereopermutator on " << e.first << " - " << e.second << ": " << stereopermutator.info() << "\n";
    }
    return false;
  }

  if(stereopermutatorOption->numStereopermutations() != numPermutations) {
    std::cout << "Bond stereopermutator on " << stereopermutatorOption->edge().first
      << " - " << stereopermutatorOption->edge().second
      << " has " << stereopermutatorOption->numStereopermutations()
      << " stereopermutations, not " << numPermutations << "\n";
    return false;
  }

  if(assignment) {
    if(stereopermutatorOption->assigned() != assignment.value()) {
      std::cout << "Bond stereopermutator on " << stereopermutatorOption->edge().first
        << " - " << stereopermutatorOption->edge().second
        << " is assigned "
        << (
          stereopermutatorOption->assigned()
          ? std::to_string(stereopermutatorOption->assigned().value())
          : "u"
        )
        << ", not " << assignment.value() << "\n";
      return false;
    }
  }

  return true;
}

bool isAtomStereocenter(
  const Molecule& molecule,
  AtomIndex i,
  const unsigned numPermutations,
  const boost::optional<unsigned>& assignment
) {
  auto stereopermutatorOption = molecule.stereopermutators().option(i);

  if(!stereopermutatorOption) {
    std::cout << "No stereopermutator on atom index " << i << "\n";
    return false;
  }

  if(stereopermutatorOption->numStereopermutations() != numPermutations) {
    std::cout << "Atom stereopermutator on " << i << " has "
      << stereopermutatorOption->numStereopermutations() << " stereopermutations, not "
      << numPermutations << "\n";
    return false;
  }

  if(assignment) {
    if(stereopermutatorOption->assigned() != assignment.value()) {
      std::cout << "Atom stereopermutator on " << i << " is assigned "
        << (
          stereopermutatorOption->assigned()
          ? std::to_string(stereopermutatorOption->assigned().value())
          : "u"
        ) << ", not " << assignment.value() << "\n";
      return false;
    }
  }

  return true;
}

bool isStereogenic(
  const Molecule& molecule,
  AtomIndex i
) {
  auto stereopermutatorOption = molecule.stereopermutators().option(i);

  if(!stereopermutatorOption) {
    return false;
  }

  if(stereopermutatorOption->numStereopermutations() <= 1) {
    return false;
  }

  return true;
}

std::string getPathString(const std::string& fileName) {
  boost::filesystem::path filePath("ranking_tree_molecules");
  filePath /= fileName;
  return filePath.string();
}

BOOST_AUTO_TEST_CASE(sequenceRuleOneTests) {
  using namespace std::string_literals;
  // P. 92.2.2 Sequence subrule 1b: Priority due to duplicate atoms
  // Cycle and multiple-bond splitting
  auto exampleThree = IO::read(
    getPathString("1S5R-bicyclo-3-1-0-hex-2-ene.mol")
  );

  auto exampleThreeExpanded = RankingTree(
    exampleThree.graph(),
    exampleThree.graph().cycles(),
    exampleThree.stereopermutators(),
    exampleThree.dumpGraphviz(),
    0,
    {},
    RankingTree::ExpansionOption::Full
  );

  auto exampleThreeRanked = exampleThreeExpanded.getRanked();

  BOOST_CHECK_MESSAGE(
    (exampleThreeRanked == std::vector<
      std::vector<unsigned long>
    > { {6}, {3}, {2}, {1} }),
    "Example three expanded is not {{6}, {3}, {2}, {1}}, but: "
    << temple::stringify(exampleThreeRanked)
  );

  auto exampleThreeExpandedAgain = RankingTree(
    exampleThree.graph(),
    exampleThree.graph().cycles(),
    exampleThree.stereopermutators(),
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

BOOST_AUTO_TEST_CASE(sequenceRuleThreeTests) {
  // P-92.4.2.1 Example 1 (Z before E)
  auto ZEDifference = IO::read(
    getPathString("2Z5S7E-nona-2,7-dien-5-ol.mol"s)
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(ZEDifference, 0, 2, 0),
    "Stereopermutator at C0 in 2Z5S7E-nona-2,7-dien-5-ol is not S"
  );

  // P-92.4.2.2 Example 1 (Z before E in aux. stereopermutators, splitting)
  auto EECyclobutane = IO::read(
    getPathString("1E3E-1,3-difluoromethylidenecyclobutane.mol")
  );

  BOOST_CHECK_MESSAGE(
    isBondStereopermutator(EECyclobutane, BondIndex {0, 3}, 2, 0)
    && isBondStereopermutator(EECyclobutane, BondIndex {5, 6}, 2, 0),
    "1E3E-1,3-difluoromethylidenecyclobutane double bonds aren't E. "
  );

  // P-92.4.2.2 Example 2 (stereogenic before non-stereogenic)
  auto inTreeNstgDB = IO::read(
    getPathString("(2Z5Z7R8Z11Z)-9-(2Z-but-2-en-1-yl)-5-(2E-but-2-en-1-yl)trideca-2,5,8,11-tetraen-7-ol.mol"s)
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(inTreeNstgDB, 0, 2, 1u),
    "(2Z5Z7R8Z11Z)-9-(2Z-but-2-en-1-yl)-5-(2E-but-2-en-1-yl)trideca-2,5,8,11-tetraen-7-ol "
    "difference between non-stereogenic auxiliary stereopermutator and assigned "
    "stereopermutator isn't recognized! "
  );
}

BOOST_AUTO_TEST_CASE(sequenceRuleFourTests) {
  // (4A) P-92.5.1 Example (stereogenic before non-stereogenic)
  auto pseudoOverNonstg = IO::read(
    getPathString("(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol"s)
  );

  BOOST_CHECK_MESSAGE(
    !isStereogenic(pseudoOverNonstg, 10),
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "branch with R-R aux. stereopermutators not non-stereogenic"
  );

  BOOST_CHECK_MESSAGE(
    isStereogenic(pseudoOverNonstg, 1),
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "branch with R-S aux. stereopermutators not stereogenic"
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(pseudoOverNonstg, 0, 2, 0),
    "(2R,3s,4S,6R)-2,6-dichloro-5-(1R-1-chloroethyl)-3-(1S-1-chloroethyl)heptan-4-ol.mol "
    "sequence rule 4A does not recognize stereogenic over non-stereogenic, 3 as S"
  );

  // (4B) P-92.5.2.2 Example 1 (single chain pairing, ordering and reference selection)
  auto simpleLikeUnlike = IO::read(
    getPathString("(2R,3R,4R,5S,6R)-2,3,4,5,6-pentachloroheptanedioic-acid.mol"s)
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(simpleLikeUnlike, 10, 2, 1u),
    "(2R,3R,4R,5S,6R)-2,3,4,5,6-pentachloroheptanedioic-acid central carbon does "
    " not register as a stereopermutator and/or isn't assigned as R"
  );

  // (4B) P-92.5.2.2 Example 3 (single-chain pairing, cycle splitting)
  auto lAlphaLindane = IO::read(
    getPathString("l-alpha-lindane.mol"s)
  );

  BOOST_CHECK_MESSAGE(
    (
      temple::all_of(
        std::vector<AtomIndex> {6, 7, 8, 9, 10, 11},
        [&](const auto carbonIndex) -> bool {
          return isStereogenic(lAlphaLindane, carbonIndex);
        }
      )
    ),
    "Not all L-alpha-lindane carbon atoms not recognized as stereopermutators!"
  );

  // (4B) P-92.5.2.2 Example 4 (multiple-chain stereopermutator ranking)
  auto oxyNitroDiffBranches = IO::read(
    getPathString("(2R,3S,6R,9R,10S)-6-chloro-5-(1R,2S)-1,2-dihydroxypropoxy-7-(1S,2S)-1,2-dihydroxypropoxy-4,8-dioxa-5,7-diazaundecande-2,3,9,10-tetrol.mol"s)
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(oxyNitroDiffBranches, 0, 2, 1u),
    "(2R,3S,6R,9R,10S)-6-chloro-5-(1R,2S)-1,2-dihydroxypropoxy-7-(1S,2S)-1,2-dihydroxypropoxy-4,8-dioxa-5,7-diazaundecande-2,3,9,10-tetrol central carbon not recognized as R"
  );

  // (4B) P-92.5.2.2 Example 5 (multiple-chain stereopermutator ranking)
  auto groupingDifferences = IO::read(
    getPathString("(2R,3R,5R,7R,8R)-4.4-bis(2S,3R-3-chlorobutan-2-yl)-6,6-bis(2S,4S-3-chlorobutan-2-yl)-2,8-dichloro-3,7-dimethylnonan-5-ol.mol"s)
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(groupingDifferences, 0, 2, 1u),
    "The central carbon in (2R,3R,5R,7R,8R)-4.4-bis(2S,3R-3-chlorobutan-2-yl)-6,6-bis(2S,4S-3-chlorobutan-2-yl)-2,8-dichloro-3,7-dimethylnonan-5-ol is not recognized as R"
  );

  // (4B) P-92.5.2.2 Example 6 (number of reference descriptors)
  auto numReferenceDescriptors = IO::read(
    getPathString("2R-2-bis(1R)-1-hydroxyethylamino-2-(1R)-1-hydroxyethyl(1S)-1-hydroxyethylaminoacetic-acid.mol")
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(numReferenceDescriptors, 0, 2, 1u),
    "The central carbon in 2R-2-bis(1R)-1-hydroxyethylamino-2-(1R)-1-hydroxyethyl(1S)-1-hydroxyethylaminoacetic-acid is not recognized as R"
  );
}

BOOST_AUTO_TEST_CASE(sequenceRuleFiveTests) {
  // (4C) P-92.5.3 Example r/s leads to R difference
  auto rsDifference = IO::read(
    getPathString("(2R,3r,4R,5s,6R)-2,6-dichloro-3,5-bis(1S-1-chloroethyl)heptan-4-ol.mol")
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(rsDifference, 0, 2, 1u),
    "The central carbon in (2R,3r,4R,5s,6R)-2,6-dichloro-3,5-bis(1S-1-chloroethyl)heptan-4-ol is not recognized as R"
  );

  // (5) P-92.6 Example 1 simple R/S difference leads to r
  auto pseudo = IO::read(
    getPathString("(2R,3r,4S)-pentane-2,3,4-trithiol.mol")
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(pseudo, 0, 2, 1u),
    "The central carbon in (2R,3r,4S)-pentane-2,3,4-trithiol is not recognized as R"
  );

  // (5) P-92.6 Example 2 cyclobutane splitting
  auto cyclobutane = IO::read(
    getPathString("(1r,3r)-cyclobutane-1,3-diol.mol")
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(cyclobutane, 2, 2, 1u)
    && isAtomStereocenter(cyclobutane, 3, 2, 1u),
    "The chiral carbons in (1r,3r)-cyclobutane-1,3-diol aren't properly recognized"
  );

  // (5) P-92.6 Example 5 double bond ranking
  auto pseudoDB = IO::read(
    getPathString("(2E,4R)-4-chloro-3-(1S-1-chloroethyl)pent-2-ene.mol")
  );

  BOOST_CHECK_MESSAGE(
    isBondStereopermutator(pseudoDB, BondIndex {0, 3}, 2, 0u),
    "Double bond in (2E,4R)-4-chloro-3-(1S-1-chloroethyl)pent-2-ene isn't E"
  );

  // (5) P-92.6 Example 6
  auto fourDoesNothing = IO::read(
    getPathString("1s-1-(1R,2R-1,2-dichloropropyl-1S,2R-1,2-dichloropropylamino)1-(1R,2S-1,2-dichloropropyl-1S,2S-1,2-dichloropropylamino)methan-1-ol.mol")
  );

  BOOST_CHECK_MESSAGE(
    isAtomStereocenter(fourDoesNothing, 0, 2, 0u),
    "The central stereopermutator in 1s-1-(1R,2R-1,2-dichloropropyl-1S,2R-1,2-dichloropropylamino)1-(1R,2S-1,2-dichloropropyl-1S,2S-1,2-dichloropropylamino)methan-1-ol isn't recognized as S"
  );
}
