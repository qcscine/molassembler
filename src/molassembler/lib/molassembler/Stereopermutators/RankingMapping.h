#ifndef INCLUDE_MOLASSEMBLER_RANKING_MAPPING_H
#define INCLUDE_MOLASSEMBLER_RANKING_MAPPING_H

#include "molassembler/RankingInformation.h"
#include <unordered_map>
#include "boost/optional/optional.hpp"

namespace Scine {
namespace molassembler {

struct SiteMapping {
  struct Changes {
    AtomIndex atom;
    SiteIndex site;
  };

  std::unordered_map<SiteIndex, SiteIndex, SiteIndex::Hash> map;
  boost::optional<Changes> changes;

  static SiteMapping from(const RankingInformation& a, const RankingInformation& b);
};

} // namespace molassembler
} // namespace Scine

#endif
