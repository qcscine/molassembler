/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Data struct storing results of ranking and local graph algorithms
 */

#ifndef INCLUDE_MOLASSEMBLER_RANKING_INFORMATION_H
#define INCLUDE_MOLASSEMBLER_RANKING_INFORMATION_H

#include "Molassembler/Types.h"
#include "Molassembler/Temple/StrongIndex.h"

#include <vector>

namespace Scine {
namespace Molassembler {

struct site_index_tag;
using SiteIndex = Temple::StrongIndex<site_index_tag, unsigned>;

//! Ranking data of substituents around a central vertex
struct MASM_EXPORT RankingInformation {
//!@name Member types
//!@{
  struct Link;

  template<typename T>
  using NestedList = std::vector<
    std::vector<T>
  >;

  //! ASC ordered list (via ranking) of atom index lists (sub-list atoms equal)
  using RankedSubstituentsType = NestedList<AtomIndex>;

  //! An unordered list of sets of atom indices that constitute binding sites
  using SiteListType = NestedList<AtomIndex>;

  //! Ascending ordered list of binding site indices (sub-list site indices equal)
  using RankedSitesType = NestedList<SiteIndex>;
//!@}

//!@name Static member functions
//!@{
  /*! @brief Gets ranking positions of a binding site's constituting atoms in
   *   descending order
   *
   * @verbatim
   *         Ranking position:     0     1       2      3
   * Input:  substituentRanking: {{0}, {4, 3}, {1, 9}, {5}}
   *         binding site indices: {5, 3, 0}
   * Output: {3, 1, 0}
   * @endverbatim
   *
   * @complexity{@math{\Theta(S)} where @math{S} is the number of substituents}
   */
  static std::vector<unsigned> siteConstitutingAtomsRankedPositions(
    const std::vector<AtomIndex>& siteAtomList,
    const RankingInformation::RankedSubstituentsType& substituentRanking
  );

  /** @brief Combines the ranking of substituents with the binding site
   *   constitutions into a site-level ranking
   *
   * Site ordering is predicated on the following rules:
   * 1. Larger sites precede smaller sites
   * 2. Lexicographical comparison of the constituting atoms' ranking positions
   *
   * @complexity{@math{\Theta(S^2)} where @math{S} is the number of sites}
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
   *   comparison operators
   */
  std::vector<Link> links;
//!@}

//!@name Modification
//!@{
  /** @brief Applies an atom index permutation
   *
   * @complexity{@math{\Theta(1)}}
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
   * @complexity{@math{\Theta(1)}}
   *
   * @throws std::out_of_range If the specified atom index is not part of any
   *   binding site
   */
  SiteIndex getSiteIndexOf(AtomIndex i) const;

  /** @brief Fetches the position of a binding site index within the site ranking
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @param i The binding site index to find
   *
   * @return The position within @p siteRanking of the supplied ligand index
   */
  unsigned getRankedIndexOfSite(SiteIndex i) const;

  /*! @brief Checks whether there are haptic binding sites
   *
   * @complexity{@math{\Theta(1)}}
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

/**
 * @brief Information on links between substituents of a central atom
 */
struct RankingInformation::Link {
//!@name Special member functions
//!@{
  /*! @brief Default constructor
   * @warning Does not establish member invariants
   */
  Link();

  /*! @brief Constructor from data without established invariants
   *
   * @complexity{@math{\Theta(N)}}
   */
  Link(
    std::pair<SiteIndex, SiteIndex> siteIndices,
    std::vector<AtomIndex> sequence,
    AtomIndex source
  );
//!@}

//!@name Data members
//!@{
  //! An (asc) ordered pair of the site indices that are linked
  std::pair<SiteIndex, SiteIndex> sites;

  /*!
   * @brief The in-order atom sequence of the cycle atom indices
   *
   * @note The cycle sequence is centralized on the source vertex, meaning the
   * front index is always the source vertex
   *
   * @note The cycle sequence is standardized by ordering the second and last
   * vertices of the sequence ascending (i.e. reversing the sequence past the
   * source vertex if the second index is larger than the last one)
   */
  std::vector<AtomIndex> cycleSequence;
//!@}

//!@name Modification
//!@{
  /*! @brief Apply an index permutation to this object. Re-establishes invariants.
   *
   * @complexity{@math{\Theta(N)}}
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);
//!@}

//!@name Operators
//!@{
  //! Performs a lexicographical comparison on both data members
  bool operator == (const Link& other) const;
  bool operator != (const Link& other) const;

  //! Performs a lexicographical comparison on both data members
  bool operator < (const Link& other) const;
//!@}
};

} // namespace Molassembler
} // namespace Scine

#endif
