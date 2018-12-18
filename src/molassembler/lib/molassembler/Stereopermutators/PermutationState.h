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
  //! Stably resorted (by set size) ligands ranking
  RankingInformation::RankedLigandsType canonicalLigands;

  //! Character representation of bonding case
  std::vector<char> symbolicCharacters;

  //! Self-referential representation of links
  stereopermutation::Stereopermutation::LinksSetType selfReferentialLinks;

  //! Mapping from ligand index to modeled ligand plane distance
  std::vector<DistanceGeometry::ValueBounds> ligandDistances;

  using ConeAngleType = std::vector<
    boost::optional<DistanceGeometry::ValueBounds>
  >;
  //! Mapping from ligand index to cone angle optional
  ConeAngleType coneAngles;

  //! Vector of rotationally unique stereopermutations with associated weights
  stereopermutation::StereopermutationsWithWeights permutations;

  //! Vector of permutation indices that are feasible
  std::vector<unsigned> feasiblePermutations;

  //! Mapping from ligand index to permutational symmetry position
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
  /*! Stably re-sort ranked ligand indices in decreasing set size
   *
   * Necessary to avoid treating e.g. AAB and ABB separately, although the
   * resulting assignments are identical.
   *
   * Example: rankedLigands: {5, 8}, {3}, {1, 2, 4}
   * Result: {1, 2, 4}, {5, 8}, {3}
   */
  static RankingInformation::RankedLigandsType canonicalize(
    RankingInformation::RankedLigandsType rankedLigands
  );

  //! Condense ligand ranking information into canonical characters for symbolic computation
  static std::vector<char> transferToSymbolicCharacters(
    const RankingInformation::RankedLigandsType& canonicalLigands
  );

  //! Make ligand-index based ligands self-referential within canonical ligands
  static stereopermutation::Stereopermutation::LinksSetType selfReferentialTransform(
    const std::vector<LinkInformation>& rankingLinks,
    const RankingInformation::RankedLigandsType& canonicalLigands
  );

  /*!
   * @brief Generates a flat mapping from ligand indices to symmetry positions
   *
   * Generates a mapping from ligand indices to symmetry positions according to
   * the ranking character distribution to symmetry positions of an assignment
   * (its characters member) and any defined links between symmetry positions.
   */
  static std::vector<unsigned> generateLigandToSymmetryPositionMap(
    const stereopermutation::Stereopermutation& assignment,
    const RankingInformation::RankedLigandsType& canonicalLigands
  );

  /*!
   * @brief Generates a flat mapping from symmetry positions to ligand indices
   *
   * Generates exactly the inverse map to generateLigandToSymmetryPositionMap
   */
  static std::vector<unsigned> generateSymmetryPositionToLigandMap(
    const stereopermutation::Stereopermutation& assignment,
    const RankingInformation::RankedLigandsType& canonicalLigands
  );

  /*!
   * @brief Generates the reduced character representation of ligands at their
   *   current symmetry positions
   */
  static std::vector<char> makeStereopermutationCharacters(
    const RankingInformation::RankedLigandsType& canonicalLigands,
    const std::vector<char>& canonicalStereopermutationCharacters,
    const std::vector<unsigned>& ligandsAtSymmetryPositions
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
    const RankingInformation::RankedLigandsType& canonicalLigands,
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
