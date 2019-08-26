/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Decide which stereopermutations are feasible
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_FEASIBLE_PERMUTATIONS_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_FEASIBLE_PERMUTATIONS_H

#include "molassembler/DistanceGeometry/ValueBounds.h"
#include "molassembler/Options.h"
#include "molassembler/RankingInformation.h"

#include "stereopermutation/Stereopermutation.h"

namespace Scine {
namespace molassembler {

// Forward-declarations
struct AbstractStereopermutations;

struct FeasibleStereopermutations {
//!@name Public types
//!@{
  using ConeAngleType = std::vector<
    boost::optional<DistanceGeometry::ValueBounds>
  >;
//!@}

//!@name Static functions
//!@{
  /*!
   * @brief Determine whether a stereopermutation is not obviously impossible
   *
   * @complexity{@math{\Theta(L)}}
   * @todo Move this to SpatialModel
   */
  static bool isNotObviouslyImpossibleStereopermutation(
    AtomIndex centralIndex,
    const stereopermutation::Stereopermutation& assignment,
    const RankingInformation::RankedSitesType& canonicalSites,
    const ConeAngleType& coneAngles,
    const RankingInformation& ranking,
    Symmetry::Name symmetry,
    const OuterGraph& graph
  );
//!@}

//!@name Constructors
//!@{
  //! Empty initializer, all data members are null objects
  FeasibleStereopermutations() = default;

  /**
   * @brief Determines the subset of stereopermutations that are feasible in
   *   three dimensions
   *
   * @param abstractPermutations The set of abstract stereopermutations
   * @param ranking Ranking object indicating chemical differences between
   *   sites and substituents
   * @param graph The graph being modeled
   *
   * @complexity{@math{\Theta(P\cdot L)} where @math{P} is the number of
   * abstract stereopermutations and @math{L} is the number of links}
   */
  FeasibleStereopermutations(
    const AbstractStereopermutations& abstractPermutations,
    Symmetry::Name symmetry,
    AtomIndex centralIndex,
    const RankingInformation& ranking,
    const OuterGraph& graph
  );
//!@}

//!@name Data members
//!@{
  //! Mapping from site index to modeled site plane distance
  std::vector<DistanceGeometry::ValueBounds> siteDistances;

  //! Mapping from site index to cone angle optional
  ConeAngleType coneAngles;

  //! Vector of permutation indices that are feasible
  std::vector<unsigned> indices;
//!@}
};

} // namespace molassembler
} // namespace Scine

#endif
