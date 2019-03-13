/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Permutational state of AtomStereopermutator member class
 */

#ifndef INCLUDE_MOLASSEMBLER_DETAIL_PERMUTATION_STATE_H
#define INCLUDE_MOLASSEMBLER_DETAIL_PERMUTATION_STATE_H

#include "stereopermutation/GenerateUniques.h"
#include "chemical_symmetries/Properties.h"

#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/Options.h"

namespace Scine {

namespace molassembler {

/**
 * @brief Contains all permutational state of an atom stereopermutator
 */
struct PermutationState {
//!@name State
//!@{
  //! Stably resorted (by set size) site ranking
  RankingInformation::RankedSitesType canonicalSites;

  //! Character representation of bonding case
  std::vector<char> symbolicCharacters;

  //! Self-referential representation of links
  stereopermutation::Stereopermutation::LinksSetType selfReferentialLinks;

  //! Mapping from site index to modeled site plane distance
  std::vector<DistanceGeometry::ValueBounds> siteDistances;

  using ConeAngleType = std::vector<
    boost::optional<DistanceGeometry::ValueBounds>
  >;
  //! Mapping from site index to cone angle optional
  ConeAngleType coneAngles;

  //! Vector of rotationally unique stereopermutations with associated weights
  stereopermutation::StereopermutationsWithWeights permutations;

  //! Vector of permutation indices that are feasible
  std::vector<unsigned> feasiblePermutations;

  //! Mapping from site index to permutational symmetry position
  std::vector<unsigned> symmetryPositionMap;
//!@}

//!@name Special member functions
//!@{
  PermutationState() = default;

  // Establish full state
  PermutationState(
    const RankingInformation& ranking,
    AtomIndex centerAtom,
    Symmetry::Name symmetry,
    const OuterGraph& graph
  );
//!@}

//!@name Static member functions
//!@{
  /*!
   * @brief Stably re-sort ranked site indices in decreasing set size
   *
   * Necessary to avoid treating e.g. AAB and ABB separately, although the
   * resulting assignments are identical.
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
   * @brief Condense site ranking information into canonical characters for
   *   symbolic computation
   *
   * Use the output of canonicalize here as input:
   *
   * Example:
   * @verbatim
   * rankedSites = {5, 8}, {3}, {1, 2, 4}
   * canonical = canonicalize(rankedSites = {1, 2, 4}, {5, 8}, {3}
   * transferToSymbolicCharacetrs(canonical) = A, A, A, B, B, C
   * @endverbatim
   */
  static std::vector<char> transferToSymbolicCharacters(
    const RankingInformation::RankedSitesType& canonicalSites
  );

  /*!
   * @brief Make site-index based links self-referential within canonical sites
   *
   * Example:
   * @verbatim
   * links = {indexPair = {5, 8}}
   *
   * self-refential idx: 0  1  2    3  4    5
   * canonicalSites = {1, 2, 4}, {5, 8}, {3} (this is output from canonicalize)
   *
   * selfReferentialTransform(links, canonicalSites) = {3, 4}
   * @endverbatim
   */
  static stereopermutation::Stereopermutation::LinksSetType selfReferentialTransform(
    const std::vector<LinkInformation>& rankingLinks,
    const RankingInformation::RankedSitesType& canonicalSites
  );

  /*!
   * @brief Generates a flat mapping from site indices to symmetry positions
   *
   * Generates a mapping from site indices to symmetry positions according to
   * the ranking character distribution to symmetry positions of an assignment
   * (its characters member) and any defined links between symmetry positions.
   *
   * @code{.cpp}
   * auto mapping = generateSiteToSymmetryPosition(...);
   * unsigned symmetryPositionOfSiteFour = mapping.at(4u);
   * @endcode
   */
  static std::vector<unsigned> generateSiteToSymmetryPositionMap(
    const stereopermutation::Stereopermutation& assignment,
    const RankingInformation::RankedSitesType& canonicalSites
  );

  /*!
   * @brief Generates a flat mapping from symmetry positions to site indices
   *
   * Generates exactly the inverse map to generateSiteToSymmetryPositionMap
   *
   * @code{cpp}
   * auto mapping = generateSymmetryPositionToSiteMap(...);
   * unsigned siteIndexAtSymmetryPositionFive = mapping.at(5u);
   * @endcode
   */
  static std::vector<unsigned> generateSymmetryPositionToSiteMap(
    const stereopermutation::Stereopermutation& assignment,
    const RankingInformation::RankedSitesType& canonicalSites
  );

  /*!
   * @brief Generates the reduced character representation of sites at their
   *   current symmetry positions
   */
  static std::vector<char> makeStereopermutationCharacters(
    const RankingInformation::RankedSitesType& canonicalSites,
    const std::vector<char>& canonicalStereopermutationCharacters,
    const std::vector<unsigned>& sitesAtSymmetryPositions
  );

  /*!
   * @brief Generates an symmetry position index mapping from a symmetry
   *   transition group
   *
   * @note This has to be a copy-initialized optional in the return value.
   *   Don't change it, unless you want to sift through -fsanitize=address
   *   output to find the bug in the optional propagation. It's not safe to
   *   return a reference to within a temporary object where this is used.
   */
  static boost::optional<
    std::vector<unsigned>
  > getIndexMapping(
    const Symmetry::properties::SymmetryTransitionGroup& mappingsGroup,
    const ChiralStatePreservation& preservationOption
  );

  /*! Determine whether a stereopermutation is not obviously impossible
   * @todo Move this to SpatialModel
   */
  static bool isNotObviouslyImpossibleStereopermutation(
    const stereopermutation::Stereopermutation& assignment,
    const RankingInformation::RankedSitesType& canonicalSites,
    const ConeAngleType& coneAngles,
    const RankingInformation& ranking,
    Symmetry::Name symmetry,
    const OuterGraph& graph
  );
//!@}
};


} // namespace molassembler

} // namespace Scine

#endif
