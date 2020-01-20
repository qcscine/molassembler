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
namespace stereopermutators {

// Forward-declarations
struct Abstract;

struct Feasible {
//!@name Public types
//!@{
  using ConeAngleType = std::vector<
    boost::optional<distance_geometry::ValueBounds>
  >;
//!@}

//!@name Static functions
//!@{
  /*! @brief Determine whether a link is possibly feasible
   *
   * Catches some obviously impossible links, but does not imply
   * that the link is truly feasible if the test passes.
   *
   * @todo Move this to SpatialModel
   */
  static bool linkPossiblyFeasible(
    const LinkInformation& link,
    AtomIndex centralIndex,
    const ConeAngleType& cones,
    const RankingInformation& ranking,
    shapes::Shape shape,
    const std::vector<unsigned>& shapeVertexMap,
    const Graph& graph
  );

  /*! @brief Determine whether a stereopermutation is possibly feasible
   *
   * Catches some obviously impossible stereopermutations, but does not
   * imply that the stereopermutation is truly feasibly if the test passes.
   *
   * @complexity{@math{\Theta(L)}}
   * @todo Move this to SpatialModel
   */
  static bool possiblyFeasible(
    const stereopermutation::Stereopermutation& assignment,
    AtomIndex centralIndex,
    const RankingInformation::RankedSitesType& canonicalSites,
    const ConeAngleType& coneAngles,
    const RankingInformation& ranking,
    shapes::Shape shape,
    const Graph& graph
  );
//!@}

//!@name Constructors
//!@{
  //! Empty initializer, all data members are null objects
  Feasible() = default;

  /**
   * @brief Determines the subset of stereopermutations that are feasible in
   *   three dimensions
   *
   * @param abstractPermutations The set of abstract stereopermutations
   * @param shape The underlying shape of the stereopermutator
   * @param centralIndex the atom index of the stereopermutator
   * @param ranking Ranking object indicating chemical differences between
   *   sites and substituents
   * @param graph The graph being modeled
   *
   * @complexity{@math{\Theta(P\cdot L)} where @math{P} is the number of
   * abstract stereopermutations and @math{L} is the number of links}
   */
  Feasible(
    const Abstract& abstractPermutations,
    shapes::Shape shape,
    AtomIndex centralIndex,
    const RankingInformation& ranking,
    const Graph& graph
  );
//!@}

//!@name Data members
//!@{
  //! Mapping from site index to modeled site plane distance
  std::vector<distance_geometry::ValueBounds> siteDistances;

  //! Mapping from site index to cone angle optional
  ConeAngleType coneAngles;

  //! Vector of permutation indices that are feasible
  std::vector<unsigned> indices;
//!@}
};

} // namespace stereopermutators
} // namespace molassembler
} // namespace Scine

#endif
