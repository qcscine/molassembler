/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "molassembler/Stereopermutators/FeasiblePermutations.h"

#include "molassembler/Stereopermutators/AbstractPermutations.h"
#include "molassembler/Stereopermutators/CycleFeasibility.h"
#include "molassembler/Stereopermutators/ShapeVertexMaps.h"
#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"

#include "shapes/PropertyCaching.h"
#include "temple/Adaptors/CyclicFrame.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"

namespace Scine {
namespace molassembler {
namespace stereopermutators {

bool Feasible::linkPossiblyFeasible(
  const LinkInformation& link,
  const AtomIndex centralIndex,
  const ConeAngleType& cones,
  const RankingInformation& ranking,
  const shapes::Shape shape,
  const SiteToShapeVertexMap& shapeVertexMap,
  const Graph& graph
) {
  // The algorithm below is explained in detail in documents/denticity_feasibility
  assert(link.cycleSequence.front() != link.cycleSequence.back());

  // Perform no checks if, for either of the sites, no cone angle could be calculated
  if(!cones.at(link.indexPair.first) || !cones.at(link.indexPair.second)) {
    return true;
  }

  const distance_geometry::ValueBounds siteIConeAngle = cones.at(link.indexPair.first).value();
  const distance_geometry::ValueBounds siteJConeAngle = cones.at(link.indexPair.second).value();

  const double symmetryAngle = distance_geometry::SpatialModel::siteCentralAngle(
    centralIndex,
    shape,
    ranking,
    shapeVertexMap,
    link.indexPair,
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
    assert(link.cycleSequence.front() == centralIndex);

    /* TODO maybe it might be better to ask in a very boolean way whether
     * the shape is willing to distort for this particular link or not
     * instead of testing the angle in such a roundabout manner...
     */

    return !stereopermutators::triangleBondTooClose(
      distance_geometry::SpatialModel::modelDistance(
        centralIndex,
        link.cycleSequence[1],
        graph.inner()
      ),
      distance_geometry::SpatialModel::modelDistance(
        centralIndex,
        link.cycleSequence[2],
        graph.inner()
      ),
      symmetryAngle,
      atom_info::bondRadius(graph.elementType(centralIndex))
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
  auto symmetryGroups = shapes::properties::positionGroups(shape);


  /* First we need to construct the cyclic polygon of the cycle sequence
   * without the central atom.
   */
  assert(link.cycleSequence.front() == centralIndex);
  auto cycleEdgeLengths = temple::map(
    temple::adaptors::cyclicFrame<2>(link.cycleSequence),
    [&](const auto& i, const auto& j) -> double {
      return distance_geometry::SpatialModel::modelDistance(i, j, graph.inner());
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

  std::vector<stereopermutators::BaseAtom> bases (1);
  auto& base = bases.front();
  base.elementType = graph.elementType(centralIndex);
  base.distanceToLeft = a;
  base.distanceToRight = b;

  auto elementTypes = temple::map(
    link.cycleSequence,
    [&](const AtomIndex i) -> Utils::ElementType {
      return graph.elementType(i);
    }
  );
  // Drop the central index's element type from this map
  elementTypes.erase(std::begin(elementTypes));

  return !stereopermutators::cycleModelContradictsGraph(
    elementTypes,
    cycleEdgeLengths,
    bases
  );
}

bool Feasible::possiblyFeasible(
  const stereopermutation::Stereopermutation& stereopermutation,
  const AtomIndex centralIndex,
  const RankingInformation::RankedSitesType& canonicalSites,
  const ConeAngleType& coneAngles,
  const RankingInformation& ranking,
  const shapes::Shape shape,
  const Graph& graph
) {
  const auto shapeVertexMap = siteToShapeVertexMap(
    stereopermutation,
    canonicalSites
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
      const double symmetryAngle = distance_geometry::SpatialModel::siteCentralAngle(
        centralIndex,
        shape,
        ranking,
        shapeVertexMap,
        {siteI, siteJ},
        graph.inner()
      );

      /* A haptic steropermutation of sites is only obviously impossible if
       * the haptic sites have no spatial freedom to arrange in a fashion
       * that does not overlap.
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
  return temple::all_of(
    ranking.links,
    [&](const auto& link) -> bool {
      return linkPossiblyFeasible(
        link,
        centralIndex,
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
  const shapes::Shape shape,
  const AtomIndex centralIndex,
  const RankingInformation& ranking,
  const Graph& graph
) {
  using ModelType = distance_geometry::SpatialModel;

  siteDistances = temple::map(
    ranking.sites,
    [&](const auto& siteAtomsList) -> distance_geometry::ValueBounds {
      return ModelType::siteDistanceFromCenter(
        siteAtomsList,
        centralIndex,
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
    || temple::sum(
      temple::adaptors::transform(
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
          centralIndex,
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
    indices = temple::iota<unsigned>(P);
  }
}

} // namespace stereopermutators
} // namespace molassembler
} // namespace Scine
