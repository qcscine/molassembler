/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/RankingInformation.h"
#include "molassembler/Temple/Stringify.h"

using namespace Scine;
using namespace Molassembler;

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

  std::vector<RankingInformation::Link> symmetricHapticPincerLinks;
  RankingInformation::Link a, b;
  a.sites = {0_s, 1_s};
  b.sites = {1_s, 2_s};

  symmetricHapticPincerLinks.push_back(a);
  symmetricHapticPincerLinks.push_back(b);

  auto symmetricHapticPincerRankedSites = RankingInformation::rankSites(
    symmetricHapticPincerLigands,
    symmetricHapticPincerRanking
  );

  BOOST_CHECK_MESSAGE(
    (symmetricHapticPincerRankedSites == RankingInformation::RankedSitesType {{0_s, 2_s}, {1_s}}),
    "Expected {{0, 2}, 1}, got " << Temple::stringifyContainer(
      symmetricHapticPincerRankedSites,
      [](const auto& siteList) -> std::string {
        return "{"s + Temple::stringifyContainer(siteList,
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

  std::vector<RankingInformation::Link> asymmetricHapticPincerLinks;
  a.sites = {0_s, 1_s};
  b.sites = {1_s, 2_s};

  asymmetricHapticPincerLinks.push_back(a);
  asymmetricHapticPincerLinks.push_back(b);

  auto asymmetricHapticPincerRankedSites = RankingInformation::rankSites(
    asymmetricHapticPincerLigands,
    asymmetricHapticPincerRanking
  );

  BOOST_CHECK((asymmetricHapticPincerRankedSites == RankingInformation::RankedSitesType {{0_s}, {2_s}, {1_s}}));
}
