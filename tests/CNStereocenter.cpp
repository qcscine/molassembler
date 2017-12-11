#define BOOST_TEST_MODULE CNStereocenterTestModule
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "CNStereocenter.h"
#include "Log.h"

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
  using namespace MoleculeManip;
  using namespace MoleculeManip::Stereocenters;
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
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        a,
        [](const auto& equalPriorityIndices) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(equalPriorityIndices) + "}"s;
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
    << TemplateMagic::condenseIterable(atomsAtPositions) << "}, got {"
    << TemplateMagic::condenseIterable(symmetryMapResult)
  );
}

BOOST_AUTO_TEST_CASE(stateCorrectness) {
  using namespace MoleculeManip;
  using namespace MoleculeManip::Stereocenters;

  Log::particulars.insert(Log::Particulars::CNStereocenterStatePropagation);

  // Create a square-pyramidal chiral center
  RankingInformation squarePyramidalRanking;
  squarePyramidalRanking.sortedSubstituents = {
    {0, 4}, {2}, {3, 5}
  };

  CNStereocenter trialStereocenter {
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
