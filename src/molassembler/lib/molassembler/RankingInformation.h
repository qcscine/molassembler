#ifndef INCLUDE_ADJACENCY_LIST_RANKING_INFORMATION_H
#define INCLUDE_ADJACENCY_LIST_RANKING_INFORMATION_H

#include "GraphAlgorithms.h"

/*! @file
 *
 * Contains the data struct declaration of ranking information.
 */

/* TODO
 * - document -> where is this used and how?
 */

namespace molassembler {

/*!
 * Data struct to store the ranking information of substituents around a
 * particular central vertex
 */
struct RankingInformation {
  using RankedType = std::vector<
    std::vector<AtomIndexType>
  >;

  using LinksType = std::vector<GraphAlgorithms::LinkInformation>;

  //! Sorted substituents grouped by priority ascending
  RankedType sortedSubstituents;

  /*!
   * Information on any links between substituents in the graph.
   * The LinkInformation struct is documented in GraphAlgorithms.h
   */
  LinksType links;
};

} // namespace molassembler

#endif
