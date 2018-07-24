#define BOOST_TEST_MODULE AtomStereocenterTestModule
#include <boost/test/unit_test.hpp>

#include "temple/Stringify.h"

#include "detail/PermutationState.h"
#include "AtomStereocenter.h"
#include "Molecule.h"
#include "Log.h"

#include <iostream>

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

BOOST_AUTO_TEST_CASE(stateCorrectness) {
  // TODO this test is no longer viable (dummyMolecule must be of appropriate
  // size and offer element types in order for isFeasiblePermutation to work)
  /*
  using namespace molassembler;
  using namespace molassembler::Stereocenters;

  Log::particulars.insert(Log::Particulars::AtomStereocenterStatePropagation);

  // Hypothetical assignment eliminations only take place if links are present, but there are no links
  Molecule dummyMolecule;

  // Create a square-pyramidal chiral center
  RankingInformation squarePyramidalRanking;
  squarePyramidalRanking.sortedSubstituents = {
    {0, 4}, {2}, {3, 5}
  };
  squarePyramidalRanking.ligands = {
    {0}, {2}, {3}, {4}, {5}
  };
  squarePyramidalRanking.ligandsRanking = RankingInformation::rankLigands(
    squarePyramidalRanking.ligands,
    squarePyramidalRanking.sortedSubstituents
  );

  AtomStereocenter trialStereocenter {
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
  octahedralRanking.ligands = {
    {0}, {1}, {2}, {3}, {4}, {5}
  };
  octahedralRanking.ligandsRanking = RankingInformation::rankLigands(
    octahedralRanking.ligands,
    octahedralRanking.sortedSubstituents
  );

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
  newSquarePyramidalRanking.ligands = {
    {0}, {1}, {2}, {3}, {4}
  };
  newSquarePyramidalRanking.ligandsRanking = RankingInformation::rankLigands(
    newSquarePyramidalRanking.ligands,
    newSquarePyramidalRanking.sortedSubstituents
  );

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
  */
}

template<typename T>
using RaggedVector = std::vector<
  std::vector<T>
>;

BOOST_AUTO_TEST_CASE(PermutationStateTests) {
  using namespace molassembler;
  using TestNamespace = molassembler::PermutationState;

  using RaggedAtoms = RaggedVector<AtomIndexType>;
  using RaggedLigands = RaggedVector<unsigned>;

  using StereopermutationPairsType = std::set<
    std::pair<unsigned, unsigned>
  >;

  auto symmetricHapticPincerRanking = RaggedAtoms {
    {1, 6},
    {2, 5},
    {3, 4}
  };

  auto symmetricHapticPincerLigands = RaggedAtoms {
    {1, 2},
    {3, 4},
    {5, 6}
  };

  std::vector<GraphAlgorithms::LinkInformation> symmetricHapticPincerLinks;
  GraphAlgorithms::LinkInformation a, b;
  a.indexPair = {0, 1};
  b.indexPair = {1, 2};

  symmetricHapticPincerLinks.push_back(std::move(a));
  symmetricHapticPincerLinks.push_back(std::move(b));

  auto symmetricHapticPincerRankedLigands = RankingInformation::rankLigands(
    symmetricHapticPincerLigands,
    symmetricHapticPincerRanking
  );

  BOOST_CHECK_MESSAGE(
    (symmetricHapticPincerRankedLigands == RaggedLigands {{0, 2}, {1}}),
    "Expected {{0, 2}, 1}, got " << temple::stringify(symmetricHapticPincerRankedLigands)
  );
  BOOST_CHECK((
    TestNamespace::transferToSymbolicCharacters(
      TestNamespace::canonicalize(symmetricHapticPincerRankedLigands)
    ) == std::vector<char> {'A', 'A', 'B'}
  ));
  BOOST_CHECK((
    TestNamespace::selfReferentialTransform(
      symmetricHapticPincerLinks,
      symmetricHapticPincerRankedLigands
    ) == StereopermutationPairsType {
      {0, 2},
      {1, 2}
    }
  ));

  auto asymmetricHapticPincerRanking = RaggedAtoms {
    {1}, {6}, {2}, {5}, {3}, {4}
  };

  auto asymmetricHapticPincerLigands = RaggedAtoms {
    {1, 2},
    {3, 4},
    {5, 6}
  };

  std::vector<GraphAlgorithms::LinkInformation> asymmetricHapticPincerLinks;
  a.indexPair = {0, 1};
  b.indexPair = {1, 2};

  asymmetricHapticPincerLinks.push_back(std::move(a));
  asymmetricHapticPincerLinks.push_back(std::move(b));

  auto asymmetricHapticPincerRankedLigands = RankingInformation::rankLigands(
    asymmetricHapticPincerLigands,
    asymmetricHapticPincerRanking
  );

  BOOST_CHECK((asymmetricHapticPincerRankedLigands == RaggedLigands {{0}, {2}, {1}}));
  BOOST_CHECK((
    TestNamespace::transferToSymbolicCharacters(
      TestNamespace::canonicalize(asymmetricHapticPincerRankedLigands)
    ) == std::vector<char> {'A', 'B', 'C'}
  ));
  BOOST_CHECK((
    TestNamespace::selfReferentialTransform(
      asymmetricHapticPincerLinks,
      asymmetricHapticPincerRankedLigands
    ) == StereopermutationPairsType {
      {0, 2},
      {1, 2}
    }
  ));
}
