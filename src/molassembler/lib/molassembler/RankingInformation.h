/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Data struct storing results of ranking and local graph algorithms
 */

#ifndef INCLUDE_MOLASSEMBLER_RANKING_INFORMATION_H
#define INCLUDE_MOLASSEMBLER_RANKING_INFORMATION_H

#include "molassembler/Types.h"

#include <vector>

namespace Scine {

namespace molassembler {

/**
 * @brief Information on links between substituents of a central atom
 */
struct LinkInformation {
//!@name Special member functions
//!@{
  /*! Default constructor
   *
   * @warning Does not establish member invariants
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
   * @note The cycle sequence is centralized on the source vertex, meaning the
   * front and back indices are the source vertex
   *
   * @note The cycle sequence between the repeated source vertices is
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

  //! ascending ordered list of ligand site indices (sub-list ligand indices equal)
  using RankedLigandsType = std::vector<
    std::vector<unsigned>
  >;
//!@}

//!@name Static member functions
//!@{
  /*!
   * @brief Gets ranking positions of a ligand's constituting atoms in descending order
   *
   * @verbatim
   *                               0     1       2      3
   * Input:  sortedSubstituents: {{0}, {4, 3}, {1, 9}, {5}}
   *         ligand: {5, 3, 0}
   * Output: {3, 1, 0}
   * @endverbatim
   */
  static std::vector<unsigned> ligandConstitutingAtomsRankedPositions(
    const std::vector<AtomIndex>& ligand,
    const RankingInformation::RankedType& sortedSubstituents
  );

  /**
   * @brief Combines the ranking of substituents with the ligand site
   *   constitutions into a ligand-level ranking
   *
   * @param ligands Ligand site constitution
   * @param sortedSubstituents The substituent-level ranking result
   *
   * @return A ranked (ascending) nested list of ligand indices
   */
  static RankedLigandsType rankLigands(
    const LigandsType& ligands,
    const RankedType& sortedSubstituents
  );
//!@}

//!@name Data members
//!@{
  //! Sorted substituents grouped by priority ascending
  RankedType sortedSubstituents;

  /*!
   * @brief Ligand atom index groupings
   *
   * This is an unordered nested list of atom index lists that each constitute
   * a ligand site. If a nested list has only a single element, the ligand is
   * 'normal'. If there are multiple elements in a nested list, the ligand site
   * is haptically bonded.
   */
  LigandsType ligands;

  /*!
   * @brief Ligands ranking
   *
   * This is a ranked (ascending) nested list of ligand indices (indices into
   * #ligands).
   */
  RankedLigandsType ligandsRanking;

  /*!
   * @brief A list of information on all links between ligand sites
   *
   * @note This list is typically sorted, enabling use of the comparator
   */
  std::vector<LinkInformation> links;
//!@}

//!@name Information
//!@{
  /*!
   * @brief Fetches the ligand index of a substituent
   *
   * @throws std::out_of_range If the specified atom index is not part of any
   *   ligand
   */
  unsigned getLigandIndexOf(AtomIndex i) const;

  /**
   * @brief Fetches the position of a ligand index within the ligands ranking
   *
   * @param i The ligand index to find
   *
   * @return The position within @p ligandsRanking of the supplied ligand index
   */
  unsigned getRankedIndexOfLigand(unsigned i) const;

  /*!
   * @brief Checks whether there are haptic ligands
   *
   * @note These are identified by the #ligands member having sets of
   *   substituents with more than one element
   */
  bool hasHapticLigands() const;
//!@}

//!@name Operators
//!@{
  bool operator == (const RankingInformation& other) const;
  bool operator != (const RankingInformation& other) const;
//!@}
};

} // namespace molassembler

} // namespace Scine

#endif
