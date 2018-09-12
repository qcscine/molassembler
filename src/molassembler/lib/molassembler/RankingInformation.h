#ifndef INCLUDE_MOLASSEMBLER_RANKING_INFORMATION_H
#define INCLUDE_MOLASSEMBLER_RANKING_INFORMATION_H

#include "molassembler/Types.h"

#include <vector>

/*! @file
 *
 * @brief Data struct storing results of ranking and local graph algorithms
 */

namespace molassembler {

struct LinkInformation {
//!@name Special member functions
//!@{
  /*! Default constructor
   *
   * \warning Does not establish member invariants
   */
  LinkInformation();

  //! Constructor from data without established invariants
  LinkInformation(
    std::pair<unsigned, unsigned> ligandIndices,
    std::vector<AtomIndex> sequence,
    AtomIndex source
  );
//!@}

//!@name Data members
//!@{
  //! An (asc) ordered pair of the ligand site indices that are linked
  std::pair<unsigned, unsigned> indexPair;

  /*! The in-order atom sequence of the cycle atom indices
   *
   * \note The cycle sequence is centralized on the source vertex, meaning the
   * front and back indices are the source vertex
   *
   * \note The cycle sequence between the repeated source vertices is
   * standardized by ordering the first and last vertices of the remaining
   * sequence ascending (i.e. reversing the part of the sequence in between
   * front and back indices if the second index is larger than the
   * second-to-last one)
   */
  std::vector<AtomIndex> cycleSequence;
//!@}

//!@name Operators
//!@{
  //! Performs a lexicographical comparison on both data members
  bool operator == (const LinkInformation& other) const;
  bool operator != (const LinkInformation& other) const;

  //! Performs a lexicographical comparison on both data members
  bool operator < (const LinkInformation& other) const;
//!@}
};

//! Ranking data of substituents around a central vertex
struct RankingInformation {
//!@name Member types
//!@{
  //! ASC ordered list (via ranking) of atom index lists (sub-list atoms equal)
  using RankedType = std::vector<
    std::vector<AtomIndex>
  >;

  //! An unordered list of sets of atom indices that constitute ligand sites
  using LigandsType = std::vector<
    std::vector<AtomIndex>
  >;

  //! ASC ordered list of ligand site indices (sub-list ligand indices equal)
  using RankedLigandsType = std::vector<
    std::vector<unsigned>
  >;

  //! A list of LinkInformation structs
  using LinksType = std::vector<LinkInformation>;
//!@}

//!@name Static member functions
//!@{
  /*! Gets ranking positions of a ligand's constituting atoms in descending order
   *
   *                               0     1       2      3
   * Input:  sortedSubstituents: {{0}, {4, 3}, {1, 9}, {5}}
   *         ligand: {5, 3, 0}
   * Output: {3, 1, 0}
   */
  static std::vector<unsigned> ligandConstitutingAtomsRankedPositions(
    const std::vector<AtomIndex>& ligand,
    const RankingInformation::RankedType& sortedSubstituents
  );

  static RankedLigandsType rankLigands(
    const LigandsType& ligands,
    const RankedType& sortedSubstituents
  );
//!@}

//!@name Data members
//!@{
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
//!@}

//!@name Information
//!@{
  unsigned getLigandIndexOf(AtomIndex i) const;

  bool hasHapticLigands() const;
//!@}

//!@name Operators
//!@{
  bool operator == (const RankingInformation& other) const;
  bool operator != (const RankingInformation& other) const;
//!@}
};

} // namespace molassembler

#endif
