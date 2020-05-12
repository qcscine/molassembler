/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Decide which stereopermutations are feasible
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_FEASIBLE_PERMUTATIONS_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_FEASIBLE_PERMUTATIONS_H

#include "Molassembler/DistanceGeometry/ValueBounds.h"
#include "Molassembler/Options.h"
#include "Molassembler/Stereopermutators/ShapeVertexMaps.h"

#include "Molassembler/Stereopermutation/Stereopermutation.h"

namespace Scine {
namespace Molassembler {
namespace Stereopermutators {

// Forward-declarations
struct Abstract;

struct Feasible {
//!@name Public types
//!@{
  using ConeAngleType = std::vector<
    boost::optional<DistanceGeometry::ValueBounds>
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
    const RankingInformation::Link& link,
    AtomIndex placement,
    const ConeAngleType& cones,
    const RankingInformation& ranking,
    Shapes::Shape shape,
    const SiteToShapeVertexMap& shapeVertexMap,
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
    const Stereopermutations::Stereopermutation& assignment,
    AtomIndex placement,
    const RankingInformation::RankedSitesType& canonicalSites,
    const ConeAngleType& coneAngles,
    const RankingInformation& ranking,
    Shapes::Shape shape,
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
   * @param placement the atom index of the stereopermutator
   * @param ranking Ranking object indicating chemical differences between
   *   sites and substituents
   * @param graph The graph being modeled
   *
   * @complexity{@math{\Theta(P\cdot L)} where @math{P} is the number of
   * abstract stereopermutations and @math{L} is the number of links}
   */
  Feasible(
    const Abstract& abstractPermutations,
    Shapes::Shape shape,
    AtomIndex placement,
    const RankingInformation& ranking,
    const Graph& graph
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

} // namespace Stereopermutators
} // namespace Molassembler
} // namespace Scine

#endif
