/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/RankingInformation.h"
#include "temple/Stringify.h"

using namespace Scine;
using namespace molassembler;

using namespace std::string_literals;

inline SiteIndex operator "" _s (unsigned long long v) { return SiteIndex(v); }

BOOST_AUTO_TEST_CASE(SiteRanking) {
  const auto symmetricHapticPincerRanking = RankingInformation::RankedSubstituentsType {
    {1, 6},
    {2, 5},
    {3, 4}
  };

  const auto symmetricHapticPincerLigands = RankingInformation::SiteListType {
    {1, 2},
    {3, 4},
    {5, 6}
  };

  std::vector<LinkInformation> symmetricHapticPincerLinks;
  LinkInformation a, b;
  a.indexPair = {0_s, 1_s};
  b.indexPair = {1_s, 2_s};

  symmetricHapticPincerLinks.push_back(a);
  symmetricHapticPincerLinks.push_back(b);

  auto symmetricHapticPincerRankedSites = RankingInformation::rankSites(
    symmetricHapticPincerLigands,
    symmetricHapticPincerRanking
  );

  BOOST_CHECK_MESSAGE(
    (symmetricHapticPincerRankedSites == RankingInformation::RankedSitesType {{0_s, 2_s}, {1_s}}),
    "Expected {{0, 2}, 1}, got " << temple::stringifyContainer(
      symmetricHapticPincerRankedSites,
      [](const auto& siteList) -> std::string {
        return "{"s + temple::stringifyContainer(siteList,
          [&](SiteIndex v) -> std::string { return std::to_string(v); }
        ) + "}"s;
      }
    )
  );

  auto asymmetricHapticPincerRanking = RankingInformation::RankedSubstituentsType {
    {1}, {6}, {2}, {5}, {3}, {4}
  };

  auto asymmetricHapticPincerLigands = RankingInformation::SiteListType {
    {1, 2},
    {3, 4},
    {5, 6}
  };

  std::vector<LinkInformation> asymmetricHapticPincerLinks;
  a.indexPair = {0_s, 1_s};
  b.indexPair = {1_s, 2_s};

  asymmetricHapticPincerLinks.push_back(a);
  asymmetricHapticPincerLinks.push_back(b);

  auto asymmetricHapticPincerRankedSites = RankingInformation::rankSites(
    asymmetricHapticPincerLigands,
    asymmetricHapticPincerRanking
  );

  BOOST_CHECK((asymmetricHapticPincerRankedSites == RankingInformation::RankedSitesType {{0_s}, {2_s}, {1_s}}));
}
