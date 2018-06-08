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
/* Typedefs */
  //! ASC ordered list (via ranking) of atom index lists (sub-list atoms equal)
  using RankedType = std::vector<
    std::vector<AtomIndexType>
  >;

  //! An unordered list of sets of atom indices that constitute ligand sites
  using LigandsType = std::vector<
    std::vector<AtomIndexType>
  >;

  //! ASC ordered list of ligand site indices (sub-list ligand indices equal)
  using RankedLigandsType = std::vector<
    std::vector<unsigned>
  >;

  //! A list of LinkInformation structs
  using LinksType = std::vector<GraphAlgorithms::LinkInformation>;


/* Static members */
  /*! Gets ranking positions of a ligand's constituting atoms in descending order
   *
   *                               0     1       2      3
   * Input:  sortedSubstituents: {{0}, {4, 3}, {1, 9}, {5}}
   *         ligand: {5, 3, 0}
   * Output: {3, 1, 0}
   */
  static std::vector<unsigned> ligandConstitutingAtomsRankedPositions(
    const std::vector<AtomIndexType>& ligand,
    const RankingInformation::RankedType& sortedSubstituents
  );

  static RankedLigandsType rankLigands(
    const LigandsType& ligands,
    const RankedType& substituentRanking
  );


/* State */
  //! Sorted substituents grouped by priority ascending
  RankedType sortedSubstituents;

  //! Ligand atom index groupings
  LigandsType ligands;

  //! Ranking of ligands
  RankedLigandsType ligandsRanking;

  /*!
   * Information on any links between substituents in the graph.
   * The LinkInformation struct is documented in GraphAlgorithms.h
   */
  LinksType links;
};

} // namespace molassembler

#endif
