/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Mapping between sites generated from two ranking information instances
 */

#ifndef INCLUDE_MOLASSEMBLER_RANKING_MAPPING_H
#define INCLUDE_MOLASSEMBLER_RANKING_MAPPING_H

#include "molassembler/RankingInformation.h"
#include <unordered_map>
#include "boost/optional/optional.hpp"

namespace Scine {
namespace Molassembler {

struct SiteMapping {
  using Map = std::unordered_map<SiteIndex, SiteIndex, SiteIndex::Hash>;

  Map map;
  boost::optional<SiteIndex> changedSite;

  static SiteMapping from(const RankingInformation& a, const RankingInformation& b);
};

} // namespace Molassembler
} // namespace Scine

#endif
