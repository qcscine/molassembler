/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Compute the set of abstract permutations
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_ABSTRACT_PERMUTATIONS_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_ABSTRACT_PERMUTATIONS_H

#include "Molassembler/Stereopermutators/ShapeVertexMaps.h"
#include "Molassembler/Stereopermutation/Manipulation.h"

namespace Scine {
namespace Molassembler {

//! @brief Stereopermutator implementation details
namespace Stereopermutators {

/**
 * @brief Class to compute the set of abstract permutations from ranking
 *   and shape
 */
struct Abstract {
//!@name Static functions
//!@{
  /*!
   * @brief Stably re-sort ranked site indices in decreasing set size
   *
   * Necessary to avoid treating e.g. AAB and ABB separately, although the
   * resulting assignments are identical.
   *
   * @complexity{@math{N \log N}}
   *
   * Example:
   * @verbatim
   * rankedSites = {5, 8}, {3}, {1, 2, 4}
   * canonicalize(rankedSites) = {1, 2, 4}, {5, 8}, {3}
   * @endverbatim
   */
  static RankingInformation::RankedSitesType canonicalize(
    RankingInformation::RankedSitesType rankedSites
  );

  /*!
   * @brief Condense site ranking information into abstract rank indices for
   *   symbolic computation
   *
   * Use the output of canonicalize here as input:
   *
   * @complexity{@math{\Theta(N)}}
   *
   * Example:
   * @verbatim
   * rankedSites = {5, 8}, {3}, {1, 2, 4}
   * canonical = canonicalize(rankedSites) = {1, 2, 4}, {5, 8}, {3}
   * transferToOccupation(canonical) = {0, 0, 0, 1, 1, 2}
   *                                (= AAABBC)
   * @endverbatim
   */
  static Stereopermutations::Stereopermutation::Occupation transferToOccupation(
    const RankingInformation::RankedSitesType& canonicalSites
  );

  /*!
   * @brief Make site-index based links shape-vertex self-referential
   *
   * @complexity{@math{\Theta(L)}}
   *
   * Example:
   * @verbatim
   * links = {sites = {5, 8}}
   *
   * self-refential idx: 0  1  2    3  4    5
   * canonicalSites =   {1, 2, 4}, {5, 8}, {3} (this is output from canonicalize)
   *
   * selfReferentialTransform(links, canonicalSites) = {3, 4}
   * @endverbatim
   */
  static Stereopermutations::Stereopermutation::OrderedLinks selfReferentialTransform(
    const std::vector<RankingInformation::Link>& rankingLinks,
    const RankingInformation::RankedSitesType& canonicalSites
  );

  /*!
   * @brief Generates the reduced character representation of sites at their
   *   current shape positions
   *
   * @complexity{@math{O(S^2)} worst case}
   */
  static Stereopermutations::Stereopermutation::Occupation makeOccupation(
    const RankingInformation::RankedSitesType& canonicalSites,
    const Stereopermutations::Stereopermutation::Occupation& canonicalOccupation,
    const Temple::StrongIndexPermutation<Shapes::Vertex, SiteIndex>& sitesAtShapeVertices
  );
//!@}

//!@name Constructors
//!@{
  //! Empty initializer, all data members are nulled
  Abstract() = default;

  /**
   * @brief Generates the set of abstract stereopermutations and intermediate
   *   data
   *
   * @complexity{The generation of permutations dominates: @math{\Theta(S!)}}
   *
   * @param ranking Ranking object indicating chemical differences between
   *    substituents and sites
   * @param shape Shape name
   */
  Abstract(
    const RankingInformation& ranking,
    Shapes::Shape shape
  );
//!@}

//!@name Data members
//!@{
  //! Stably resorted (by set size) site ranking
  RankingInformation::RankedSitesType canonicalSites;

  //! Abstract representation of bonding case
  Stereopermutations::Stereopermutation::Occupation occupation;

  //! Self-referential representation of links
  Stereopermutations::Stereopermutation::OrderedLinks selfReferentialLinks;

  //! Vector of rotationally unique stereopermutations with associated weights
  Stereopermutations::Uniques permutations;
//!@}
};

} // namespace Stereopermutators
} // namespace Molassembler
} // namespace Scine


#endif
