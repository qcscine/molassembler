#ifndef INCLUDE_ADJACENCY_LIST_RANKING_INFORMATION_H
#define INCLUDE_ADJACENCY_LIST_RANKING_INFORMATION_H

#include <vector>
#include <set>
#include "common_typedefs.h"

/*! @file
 *
 * Contains the data struct declaration of ranking information.
 */

/* TODO
 * - document -> where is this used and how?
 */

namespace MoleculeManip {

/*!
 * Data struct to store the ranking information of substituents around a
 * particular central vertex
 */
struct RankingInformation {
  using RankedType = std::vector<
    std::vector<AtomIndexType>
  >;

  using LinksType = std::set<
    std::pair<AtomIndexType, AtomIndexType>
  >;

  //! Sorted substituents grouped by priority ascending
  RankedType sortedSubstituents;

  /*!
   * The set containing pairs of substituents that are linked and therefore form
   * multidentate ligands.
   */
  LinksType linkedPairs;
};

} // namespace MoleculeManip

#endif
