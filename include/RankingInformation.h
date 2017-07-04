#ifndef INCLUDE_ADJACENCY_LIST_RANKING_INFORMATION_H
#define INCLUDE_ADJACENCY_LIST_RANKING_INFORMATION_H

#include <vector>
#include <set>
#include "common_typedefs.h"

namespace MoleculeManip {

struct RankingInformation {
  std::vector<AtomIndexType> sortedPriorities;
  std::set<
    std::pair<AtomIndexType, AtomIndexType>
  > equalPriorityPairsSet;
  std::set<
    std::pair<AtomIndexType, AtomIndexType>
  > linkedPairsSet;
};

} // namespace MoleculeManip

#endif
