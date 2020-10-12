/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "Molassembler/Stereopermutators/FeasiblePermutations.h"

#include "Molassembler/Stereopermutators/AbstractPermutations.h"
#include "Molassembler/Stereopermutators/CycleFeasibility.h"
#include "Molassembler/Stereopermutators/ShapeVertexMaps.h"
#include "Molassembler/Modeling/BondDistance.h"
#include "Molassembler/Modeling/CommonTrig.h"
#include "Molassembler/DistanceGeometry/SpatialModel.h"

#include "Molassembler/Shapes/PropertyCaching.h"
#include "Molassembler/Temple/Adaptors/CyclicFrame.h"
#include "Molassembler/Temple/Adaptors/Transform.h"
#include "Molassembler/Temple/Functional.h"

namespace Scine {
namespace Molassembler {
namespace Stereopermutators {

bool Feasible::linkPossiblyFeasible(
  const RankingInformation::Link& link,
  const AtomIndex placement,
  const ConeAngleType& cones,
  const RankingInformation& ranking,
  const Shapes::Shape shape,
  const SiteToShapeVertexMap& shapeVertexMap,
  const Graph& graph
) {
  // The algorithm below is explained in detail in documents/denticity_feasibility
  assert(link.cycleSequence.front() != link.cycleSequence.back());

  // Perform no checks if, for either of the sites, no cone angle could be calculated
  if(!cones.at(link.sites.first) || !cones.at(link.sites.second)) {
    return true;
  }

  const DistanceGeometry::ValueBounds siteIConeAngle = cones.at(link.sites.first).value();
  const DistanceGeometry::ValueBounds siteJConeAngle = cones.at(link.sites.second).value();

  const double symmetryAngle = DistanceGeometry::SpatialModel::siteCentralAngle(
    placement,
    shape,
    ranking,
    shapeVertexMap,
    link.sites,
    graph.inner()
  );

  if(link.cycleSequence.size() == 3) {
    /* Test if the bond gets closer than the central index's bond radius.
     * It's a pretty arbitrary condition...
     *
     * This way, triangles aren't universally allowed (siteCentralAngle will
     * distort the angle to enable the graph in some situations, and leave
     * the ideal angle preserved in others)
     */
    assert(link.cycleSequence.front() == placement);

    /* TODO maybe it might be better to ask in a very boolean way whether
     * the shape is willing to distort for this particular link or not
     * instead of testing the angle in such a roundabout manner...
     */

    return !Stereopermutators::triangleBondTooClose(
      DistanceGeometry::SpatialModel::modelDistance(
        placement,
        link.cycleSequence[1],
        graph.inner()
      ),
      DistanceGeometry::SpatialModel::modelDistance(
        placement,
        link.cycleSequence[2],
        graph.inner()
      ),
      symmetryAngle,
      AtomInfo::bondRadius(graph.elementType(placement))
    );
  }

  /* A link across haptic sites is only obviously impossible if it is
   * impossible in the best case scenario. In this case, especially for alpha,
   * site bridge links must be possible only in the best case spatial
   * arrangement for the haptic site link to be possible. That means
   * subtracting the upper bound of the respective cone angles.
   */
  const double alpha = std::max(
    0.0,
    symmetryAngle - siteIConeAngle.upper - siteJConeAngle.upper
  );

  /* We need to respect the graph as ground truth. If a cycle is of size
   * three, then angles of whatever shape is present will be distorted if
   * there is only one group of shape positions (e.g. tetrahedral has
   * only one group, but square pyramidal has two!) available.
   *
   * For symmetries with multiple groups of shape positions, only that
   * group of shape positions with lower cross-angles is viable for the
   * link, and will distort accordingly.
   */
  // auto symmetryGroups = Shapes::Properties::positionGroupCharacters(shape);


  /* First we need to construct the cyclic polygon of the cycle sequence
   * without the central atom.
   */
  assert(link.cycleSequence.front() == placement);
  auto cycleEdgeLengths = Temple::map(
    Temple::Adaptors::cyclicFrame<2>(link.cycleSequence),
    [&](const auto& i, const auto& j) -> double {
      return DistanceGeometry::SpatialModel::modelDistance(i, j, graph.inner());
    }
  );

  /* The first and last cycle edge lengths are from and to the central atom,
   * and so we remove those by combining alpha with those edge lengths
   * with the law of cosines
   */
  const double a = cycleEdgeLengths.front();
  const double b = cycleEdgeLengths.back();
  const double c = CommonTrig::lawOfCosines(a, b, alpha); // B-A

  cycleEdgeLengths.back() = c;
  cycleEdgeLengths.erase(std::begin(cycleEdgeLengths));

  std::vector<Stereopermutators::BaseAtom> bases (1);
  auto& base = bases.front();
  base.elementType = graph.elementType(placement);
  base.distanceToLeft = a;
  base.distanceToRight = b;

  auto elementTypes = Temple::map(
    link.cycleSequence,
    [&](const AtomIndex i) -> Utils::ElementType {
      return graph.elementType(i);
    }
  );
  // Drop the central index's element type from this map
  elementTypes.erase(std::begin(elementTypes));

  return !Stereopermutators::cycleModelContradictsGraph(
    elementTypes,
    cycleEdgeLengths,
    bases
  );
}

bool Feasible::possiblyFeasible(
  const Stereopermutations::Stereopermutation& stereopermutation,
  const AtomIndex placement,
  const RankingInformation::RankedSitesType& canonicalSites,
  const ConeAngleType& coneAngles,
  const RankingInformation& ranking,
  const Shapes::Shape shape,
  const Graph& graph
) {
  const auto shapeVertexMap = siteToShapeVertexMap(
    stereopermutation,
    canonicalSites,
    ranking.links
  );

  // Check if any haptic site cones intersect
  const unsigned L = ranking.sites.size();
  for(SiteIndex siteI {0}; siteI < L - 1; ++siteI) {
    if(ranking.sites.at(siteI).size() == 1) {
      continue;
    }

    for(SiteIndex siteJ {siteI + 1}; siteJ < L; ++siteJ) {
      if(ranking.sites.at(siteJ).size() == 1) {
        continue;
      }

      // Do not test cone angles if no angle could be calculated
      if(!coneAngles.at(siteI) || !coneAngles.at(siteJ)) {
        continue;
      }

      // siteCentralAngle yields undistorted symmetry angles for haptic sites
      const double symmetryAngle = DistanceGeometry::SpatialModel::siteCentralAngle(
        placement,
        shape,
        ranking,
        shapeVertexMap,
        {siteI, siteJ},
        graph.inner()
      );

      /* A haptic steropermutation of sites is only feasible if the haptic
       * sites have spatial freedom to arrange in a fashion that does not
       * overlap.
       */
      if(
        (
          symmetryAngle
          - coneAngles.at(siteI).value().lower
          - coneAngles.at(siteJ).value().lower
        ) < 0
      ) {
        return false;
      }
    }
  }

  /* Idea: An stereopermutation is possibly feasible if all links' cycles can
   * be realized as a flat cyclic polygon, in which the edges from the central
   * atom are merged using the joint angle calculable from the
   * stereopermutation and shape.
   */
  return Temple::all_of(
    ranking.links,
    [&](const auto& link) -> bool {
      return linkPossiblyFeasible(
        link,
        placement,
        coneAngles,
        ranking,
        shape,
        shapeVertexMap,
        graph
      );
    }
  );
}

Feasible::Feasible(
  const Abstract& abstractPermutations,
  const Shapes::Shape shape,
  const AtomIndex placement,
  const RankingInformation& ranking,
  const Graph& graph
) {
  using ModelType = DistanceGeometry::SpatialModel;

  siteDistances = Temple::map(
    ranking.sites,
    [&](const auto& siteAtomsList) -> DistanceGeometry::ValueBounds {
      return ModelType::siteDistanceFromCenter(
        siteAtomsList,
        placement,
        graph
      );
    }
  );

  coneAngles.reserve(ranking.sites.size());

  for(unsigned i = 0; i < ranking.sites.size(); ++i) {
    coneAngles.push_back(
      ModelType::coneAngle(
        ranking.sites.at(i),
        siteDistances.at(i),
        graph
      )
    );
  }

  // Determine which permutations are feasible and which aren't
  const unsigned P = abstractPermutations.permutations.list.size();
  if(
    // Links are present
    !ranking.links.empty()
    // OR there are haptic sites
    || Temple::sum(
      Temple::Adaptors::transform(
        ranking.sites,
        [](const auto& siteAtomsList) -> unsigned {
          if(siteAtomsList.size() == 1) {
            return 0;
          }

          return 1;
        }
      )
    ) > 0
  ) {
    indices.reserve(P);
    for(unsigned i = 0; i < P; ++i) {
      if(
        possiblyFeasible(
          abstractPermutations.permutations.list.at(i),
          placement,
          abstractPermutations.canonicalSites,
          coneAngles,
          ranking,
          shape,
          graph
        )
      ) {
        indices.push_back(i);
      }
    }
    indices.shrink_to_fit();
  } else {
    indices = Temple::iota<unsigned>(P);
  }
}

} // namespace Stereopermutators
} // namespace Molassembler
} // namespace Scine
