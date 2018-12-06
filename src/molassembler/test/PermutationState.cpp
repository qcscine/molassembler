/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Stereopermutators/PermutationState.h"

#include "temple/Stringify.h"

template<typename T>
using RaggedVector = std::vector<
  std::vector<T>
>;

BOOST_AUTO_TEST_CASE(PermutationStateTests) {
  using namespace molassembler;
  using TestNamespace = molassembler::PermutationState;

  using RaggedAtoms = RaggedVector<AtomIndex>;
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

  std::vector<LinkInformation> symmetricHapticPincerLinks;
  LinkInformation a, b;
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

  std::vector<LinkInformation> asymmetricHapticPincerLinks;
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
