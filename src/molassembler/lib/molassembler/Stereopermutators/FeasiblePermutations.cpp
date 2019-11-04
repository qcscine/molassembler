/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
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

bool FeasibleStereopermutations::isNotObviouslyImpossibleStereopermutation(
  const AtomIndex centralIndex,
  const stereopermutation::Stereopermutation& stereopermutation,
  const RankingInformation::RankedSitesType& canonicalSites,
  const ConeAngleType& coneAngles,
  const RankingInformation& ranking,
  const Symmetry::Shape shape,
  const OuterGraph& graph
) {
  auto shapeVertexMap = siteToShapeVertexMap(
    stereopermutation,
    canonicalSites
  );

  // Check if any haptic site cones intersect
  const unsigned L = ranking.sites.size();
  for(unsigned siteI = 0; siteI < L - 1; ++siteI) {
    if(ranking.sites[siteI].size() == 1) {
      continue;
    }

    for(unsigned siteJ = siteI + 1; siteJ < L; ++siteJ) {
      if(ranking.sites[siteJ].size() == 1) {
        continue;
      }

      // Do not test cone angles if no angle could be calculated
      if(!coneAngles.at(siteI) || !coneAngles.at(siteJ)) {
        continue;
      }

      // siteCentralAngle yields undistorted symmetry angles for haptic sites
      double symmetryAngle = DistanceGeometry::SpatialModel::siteCentralAngle(
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

  /* Idea: An stereopermutation is unfeasible if any link's cycle cannot be realized
   * as a flat cyclic polygon, in which the edges from the central atom are
   * merged using the joint angle calculable from the stereopermutation and shape.
   *
   * The algorithm below is explained in detail in
   * documents/denticity_feasibility/.
   */
  for(const auto& link : ranking.links) {
    // Ignore cycles of size 3
    assert(link.cycleSequence.front() != link.cycleSequence.back());

    // Perform no checks if, for either of the sites, no cone angle could be calculated
    if(!coneAngles.at(link.indexPair.first) || !coneAngles.at(link.indexPair.second)) {
      continue;
    }

    const DistanceGeometry::ValueBounds siteIConeAngle = coneAngles.at(link.indexPair.first).value();
    const DistanceGeometry::ValueBounds siteJConeAngle = coneAngles.at(link.indexPair.second).value();

    const double symmetryAngle = DistanceGeometry::SpatialModel::siteCentralAngle(
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

      return !Stereopermutators::triangleBondTooClose(
        DistanceGeometry::SpatialModel::modelDistance(
          centralIndex,
          link.cycleSequence[1],
          graph.inner()
        ),
        DistanceGeometry::SpatialModel::modelDistance(
          centralIndex,
          link.cycleSequence[2],
          graph.inner()
        ),
        symmetryAngle,
        AtomInfo::bondRadius(graph.elementType(centralIndex))
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
    auto symmetryGroups = Symmetry::properties::positionGroups(shape);


    /* First we need to construct the cyclic polygon of the cycle sequence
     * without the central atom.
     */
    assert(link.cycleSequence.front() == centralIndex);
    auto cycleEdgeLengths = temple::map(
      temple::adaptors::cyclicFrame<2>(link.cycleSequence),
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
    cycleEdgeLengths.erase(cycleEdgeLengths.begin());

    Stereopermutators::BaseAtom base {
      graph.elementType(centralIndex),
      a,
      b,
      0.0
    };

    auto elementTypes = temple::map(
      link.cycleSequence,
      [&](const AtomIndex i) -> Utils::ElementType {
        return graph.elementType(i);
      }
    );
    // Drop the central index's element type from this map
    elementTypes.erase(elementTypes.begin());

    if(Stereopermutators::cycleModelContradictsGraph(
      elementTypes,
      cycleEdgeLengths,
      {base}
    )) {
      return false;
    }
  }

  return true;
}

FeasibleStereopermutations::FeasibleStereopermutations(
  const AbstractStereopermutations& abstractPermutations,
  const Symmetry::Shape shape,
  const AtomIndex centralIndex,
  const RankingInformation& ranking,
  const OuterGraph& graph
) {
  using ModelType = DistanceGeometry::SpatialModel;

  siteDistances = temple::map(
    ranking.sites,
    [&](const auto& siteAtomsList) -> DistanceGeometry::ValueBounds {
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
  const unsigned P = abstractPermutations.permutations.stereopermutations.size();
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
        isNotObviouslyImpossibleStereopermutation(
          centralIndex,
          abstractPermutations.permutations.stereopermutations.at(i),
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

} // namespace molassembler
} // namespace Scine
