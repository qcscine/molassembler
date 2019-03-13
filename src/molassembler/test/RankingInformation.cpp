/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/RankingInformation.h"
#include "temple/Stringify.h"

template<typename T>
using RaggedVector = std::vector<
  std::vector<T>
>;

BOOST_AUTO_TEST_CASE(rankingCombinationTests) {
  using namespace Scine;
  using namespace molassembler;

  using RaggedAtoms = RaggedVector<AtomIndex>;
  using RaggedLigands = RaggedVector<unsigned>;

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

  symmetricHapticPincerLinks.push_back(a);
  symmetricHapticPincerLinks.push_back(b);

  auto symmetricHapticPincerRankedLigands = RankingInformation::rankSites(
    symmetricHapticPincerLigands,
    symmetricHapticPincerRanking
  );

  BOOST_CHECK_MESSAGE(
    (symmetricHapticPincerRankedLigands == RaggedLigands {{0, 2}, {1}}),
    "Expected {{0, 2}, 1}, got " << temple::stringify(symmetricHapticPincerRankedLigands)
  );

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

  asymmetricHapticPincerLinks.push_back(a);
  asymmetricHapticPincerLinks.push_back(b);

  auto asymmetricHapticPincerRankedLigands = RankingInformation::rankSites(
    asymmetricHapticPincerLigands,
    asymmetricHapticPincerRanking
  );

  BOOST_CHECK((asymmetricHapticPincerRankedLigands == RaggedLigands {{0}, {2}, {1}}));
}
