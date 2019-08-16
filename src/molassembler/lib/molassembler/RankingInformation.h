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
  /*!
   * @brief Default constructor
   * @warning Does not establish member invariants
   */
  LinkInformation();

  //! Constructor from data without established invariants
  LinkInformation(
    std::pair<unsigned, unsigned> siteIndices,
    std::vector<AtomIndex> sequence,
    AtomIndex source
  );
//!@}

//!@name Data members
//!@{
  //! An (asc) ordered pair of the site indices that are linked
  std::pair<unsigned, unsigned> indexPair;

  /*!
   * @brief The in-order atom sequence of the cycle atom indices
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

//!@name Modification
//!@{
  //! Apply an index permutation to this object. Re-establishes invariants.
  void applyPermutation(const std::vector<AtomIndex>& permutation);
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
  template<typename T>
  using NestedList = std::vector<
    std::vector<T>
  >;

  //! ASC ordered list (via ranking) of atom index lists (sub-list atoms equal)
  using RankedSubstituentsType = NestedList<AtomIndex>;

  //! An unordered list of sets of atom indices that constitute binding sites
  using SiteListType = NestedList<AtomIndex>;

  //! Ascending ordered list of binding site indices (sub-list site indices equal)
  using RankedSitesType = NestedList<unsigned>;
//!@}

//!@name Static member functions
//!@{
  /*!
   * @brief Gets ranking positions of a binding site's constituting atoms in
   *   descending order
   *
   * @verbatim
   *                               0     1       2      3
   * Input:  substituentRanking: {{0}, {4, 3}, {1, 9}, {5}}
   *         binding site indices: {5, 3, 0}
   * Output: {3, 1, 0}
   * @endverbatim
   */
  static std::vector<unsigned> siteConstitutingAtomsRankedPositions(
    const std::vector<AtomIndex>& siteAtomList,
    const RankingInformation::RankedSubstituentsType& substituentRanking
  );

  /**
   * @brief Combines the ranking of substituents with the binding site
   *   constitutions into a site-level ranking
   *
   * Site ordering is predicated on the following rules:
   * 1. Larger sites precede smaller sites
   * 2. Lexicographical comparison of the constituting atoms' ranking positions
   *
   * @param sites List of atom index lists that each constitute a binding site
   * @param substituentRanking The substituent-level ranking result
   *
   * @return A ranked (ascending) nested list of ligand indices
   */
  static RankedSitesType rankSites(
    const SiteListType& sites,
    const RankedSubstituentsType& substituentRanking
  );
//!@}

//!@name Data members
//!@{
  //! Sorted substituents grouped by priority ascending
  RankedSubstituentsType substituentRanking;

  /*!
   * @brief Binding site atom index groupings
   *
   * This is an unordered nested list of atom index lists that each constitute
   * a binding site. If a nested list has only a single element, the site is
   * 'normal'. If there are multiple elements in a nested list, the site
   * is haptically bonded.
   */
  SiteListType sites;

  /*!
   * @brief Site ranking
   *
   * This is a ranked (ascending) nested list of site indices (indices into
   * #sites).
   *
   * For any instance @p r, this member should merely be the cached result of a
   * call to @p RankingInformation::rankSites with the arguments @p r.sites and
   * @p r.substituentRanking.
   */
  RankedSitesType siteRanking;

  /*!
   * @brief A list of information on all links between binding sites
   *
   * @note This list is sorted, enabling use of the RankingInformation
   *   comparatison operators
   */
  std::vector<LinkInformation> links;
//!@}

//!@name Modification
//!@{
  /**
   * @brief Applies an atom index permutation
   *
   * @param permutation The permutation to apply
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);
//!@}

//!@name Information
//!@{
  /*!
   * @brief Fetches the binding site index of a substituent
   *
   * @throws std::out_of_range If the specified atom index is not part of any
   *   binding site
   */
  unsigned getSiteIndexOf(AtomIndex i) const;

  /**
   * @brief Fetches the position of a binding site index within the site ranking
   *
   * @param i The binding site index to find
   *
   * @return The position within @p siteRanking of the supplied ligand index
   */
  unsigned getRankedIndexOfSite(unsigned i) const;

  /*!
   * @brief Checks whether there are haptic binding sites
   *
   * @note These are identified by the #sites member having sets of
   *   substituents with more than one element
   */
  bool hasHapticSites() const;
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
