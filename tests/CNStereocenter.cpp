#define BOOST_TEST_MODULE CNStereocenterTestModule
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "Molecule.h"
#include "CNStereocenter.h"
#include "Log.h"
#include "temple/Stringify.h"

std::string makeString(
  const std::vector<char>& charVec
) {
  std::string a;

  for(const auto& character : charVec) {
    a += character;
  }

  return a;
}

template<typename T, typename U>
std::string condenseMap(const std::map<T, U>& map) {
  using namespace std::string_literals;

  auto lastElementIt = map.end();
  --lastElementIt;

  std::string condensed;
  for(auto it = map.begin(); it != map.end(); ++it) {
    condensed += "{"s 
      + std::to_string(it->first) 
      + " -> "s
      + std::to_string(it->second)
      + "}"s;

    if(it != lastElementIt) {
      condensed += ", "s;
    }

  }
  return condensed;
}

BOOST_AUTO_TEST_CASE(glueTests) {
  using namespace molassembler;
  using namespace molassembler::Stereocenters;
  using namespace std::string_literals;

  // Canonicalization and canon-char generation
  // AAB vs ABB
  RankingInformation::RankedType a {{0}, {1, 2}, {3}, {4, 5}};
  auto expectedChars = glue::makeCanonicalCharacters(
    glue::canonicalize(a)
  );

  bool pass = true;
  do {
    if(glue::makeCanonicalCharacters(glue::canonicalize(a)) != expectedChars) {
      pass = false;
      break;
    }
  } while(std::next_permutation(a.begin(), a.end()));

  BOOST_CHECK_MESSAGE(
    pass,
    "Combination of canonicalization and character generation is irregular. "
    "Got an unexpected character sequence for the permutation "
    << temple::condenseIterable(
      temple::map(
        a,
        [](const auto& equalPriorityIndices) -> std::string {
          return "{"s + temple::condenseIterable(equalPriorityIndices) + "}"s;
        }
      )
    ) << " -> " << makeString(glue::makeCanonicalCharacters(glue::canonicalize(a)))
  );

  // Symmetry position mapping
  RankingInformation::RankedType rankedVariety {
    {0, 4}, {2}, {3, 5}, {1}
  };
  std::vector<AtomIndexType> atomsAtPositions {3, 1, 5, 0, 4, 2};

  auto canonRanked = glue::canonicalize(rankedVariety);
  BOOST_CHECK((canonRanked == RankingInformation::RankedType {{0, 4}, {3, 5}, {2}, {1}}));

  auto canonCharacters = glue::makeCanonicalCharacters(canonRanked);
  BOOST_CHECK((canonCharacters == std::vector<char> {'A', 'A', 'B', 'B', 'C', 'D'}));

  auto assignmentCharacters = glue::makeAssignmentCharacters(
    canonRanked,
    canonCharacters,
    atomsAtPositions
  );
  BOOST_CHECK((assignmentCharacters == std::vector<char> {'B', 'D', 'B', 'A', 'A', 'C'}));

  UniqueAssignments::Assignment sampleOctahedral {
    Symmetry::Name::Octahedral,
    assignmentCharacters
  };

  auto symmetryPositionMap = glue::makeSymmetryPositionMap(sampleOctahedral, canonRanked);
  std::map<AtomIndexType, unsigned> expectedMap {
    {0, 3},
    {1, 1},
    {2, 5},
    {3, 0},
    {4, 4},
    {5, 2}
  };
  BOOST_CHECK_MESSAGE(
    symmetryPositionMap == expectedMap,
    "makeSymmetryPositionmap returns an unexpected result! Expected {"
    << condenseMap(expectedMap) << "}, got {" << condenseMap(symmetryPositionMap)
  );

  auto symmetryMapResult = glue::mapToSymmetryPositions(sampleOctahedral, canonRanked);
  BOOST_CHECK_MESSAGE(
    symmetryMapResult == atomsAtPositions,
    "mapToSymmetryPositions returns an unexpected result! Expected {"
    << temple::condenseIterable(atomsAtPositions) << "}, got {"
    << temple::condenseIterable(symmetryMapResult)
  );

  // makeSymmetryPositionMap isn't good enough yet
  UniqueAssignments::Assignment b {
    Symmetry::Name::Octahedral,
    {'A', 'A', 'A', 'A', 'A', 'A'},
    {
      {0, 1},
      {2, 3},
      {4, 5}
    }
  };

  UniqueAssignments::Assignment c {
    Symmetry::Name::Octahedral,
    {'A', 'A', 'A', 'A', 'A', 'A'},
    {
      {0, 1},
      {2, 4},
      {3, 5}
    }
  };

  RankingInformation::RankedType canonBC {{10, 11, 12, 13, 14, 15}};

  auto bMap = glue::makeSymmetryPositionMap(b, canonBC);
  auto cMap = glue::makeSymmetryPositionMap(c, canonBC);

  BOOST_CHECK(bMap != cMap);
}

BOOST_AUTO_TEST_CASE(stateCorrectness) {
  using namespace molassembler;
  using namespace molassembler::Stereocenters;

  Log::particulars.insert(Log::Particulars::CNStereocenterStatePropagation);

  // Hypothetical assignment eliminations only take place if links are present, but there are no links
  Molecule dummyMolecule;

  // Create a square-pyramidal chiral center
  RankingInformation squarePyramidalRanking;
  squarePyramidalRanking.sortedSubstituents = {
    {0, 4}, {2}, {3, 5}
  };

  CNStereocenter trialStereocenter {
    dummyMolecule,
    Symmetry::Name::SquarePyramidal,
    8,
    squarePyramidalRanking
  };

  trialStereocenter.assign(0u);

  // Add a substituent up to octahedral
  RankingInformation octahedralRanking;
  octahedralRanking.sortedSubstituents = {
    {0, 4}, {2}, {3, 5}, {1}
  };

  trialStereocenter.addSubstituent(
    dummyMolecule,
    1,
    octahedralRanking,
    Symmetry::Name::Octahedral,
    ChiralStatePreservation::EffortlessAndUnique
  );

  BOOST_CHECK_MESSAGE(
    trialStereocenter.assigned() != boost::none,
    "Square pyramidal to Octahedral substituent addition does not preserve chiral information!"
  );

  // Simulate that the added substituent gets deleted
  trialStereocenter.propagateVertexRemoval(1);

  // And now notified that it has lost a substituent
  RankingInformation newSquarePyramidalRanking;
  newSquarePyramidalRanking.sortedSubstituents = {
    {0, 3}, {1}, {2, 4}
  };

  trialStereocenter.removeSubstituent(
    dummyMolecule,
    std::numeric_limits<AtomIndexType>::max(),
    newSquarePyramidalRanking,
    Symmetry::Name::SquarePyramidal,
    ChiralStatePreservation::EffortlessAndUnique
  );

  BOOST_CHECK_MESSAGE(
    trialStereocenter.assigned() != boost::none,
    "Octahedral to Square-pyramidal substituent removal does not preserve chiral information!"
  );

  BOOST_CHECK_MESSAGE(
    trialStereocenter.assigned() == 0u,
    "Addition and removal consistency check fails: Initial assignment is not recovered!"
  );
}

BOOST_AUTO_TEST_CASE(adhesiveTests) {
  using namespace molassembler::Stereocenters::adhesive;
  using namespace molassembler;
  using NestedVector = std::vector<
    std::vector<AtomIndexType>
  >;
  using AssignmentPairsType = std::set<
    std::pair<unsigned, unsigned>
  >;

  auto symmetricHapticPincerRanking = NestedVector {
    {1, 6},
    {2, 5},
    {3, 4}
  };

  auto symmetricHapticPincerLigands = NestedVector {
    {1, 2},
    {3, 4},
    {5, 6}
  };

  std::vector<GraphAlgorithms::LinkInformation> symmetricHapticPincerLinks;
  GraphAlgorithms::LinkInformation a, b;
  a.indexPair = {2, 3};
  b.indexPair = {4, 5};

  symmetricHapticPincerLinks.push_back(std::move(a));
  symmetricHapticPincerLinks.push_back(std::move(b));

  auto symmetricHapticPincerRankedLigands = ligandRanking(
    symmetricHapticPincerRanking,
    symmetricHapticPincerLigands
  );

  BOOST_CHECK((symmetricHapticPincerRankedLigands == NestedVector {{0, 2}, {1}}));
  BOOST_CHECK((
    canonicalCharacters(symmetricHapticPincerRankedLigands) 
    == std::vector<char> {'A', 'A', 'B'}
    ));
  BOOST_CHECK((
    canonicalLinks(
      symmetricHapticPincerLigands,
      symmetricHapticPincerRankedLigands,
      symmetricHapticPincerLinks
    ) == AssignmentPairsType {
      {0, 2},
      {1, 2}
    }
  ));

  auto asymmetricHapticPincerRanking = NestedVector {
    {1}, {6}, {2}, {5}, {3}, {4}
  };

  auto asymmetricHapticPincerLigands = NestedVector {
    {1, 2},
    {3, 4},
    {5, 6}
  };

  std::vector<GraphAlgorithms::LinkInformation> asymmetricHapticPincerLinks;
  a.indexPair = {2, 3};
  b.indexPair = {4, 5};

  asymmetricHapticPincerLinks.push_back(std::move(a));
  asymmetricHapticPincerLinks.push_back(std::move(b));

  auto asymmetricHapticPincerRankedLigands = ligandRanking(
    asymmetricHapticPincerRanking,
    asymmetricHapticPincerLigands
  );

  BOOST_CHECK((asymmetricHapticPincerRankedLigands == NestedVector {{0}, {2}, {1}}));
  BOOST_CHECK((
    canonicalCharacters(asymmetricHapticPincerRankedLigands) 
    == std::vector<char> {'A', 'B', 'C'}
    ));
  BOOST_CHECK((
    canonicalLinks(
      asymmetricHapticPincerLigands,
      asymmetricHapticPincerRankedLigands,
      asymmetricHapticPincerLinks
    ) == AssignmentPairsType {
      {0, 2},
      {1, 2}
    }
  ));
}
