/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Stereopermutators/BondStereopermutatorImpl.h"

#include "molassembler/Shapes/Data.h"

#include "molassembler/Temple/Adaptors/AllPairs.h"
#include "molassembler/Temple/Adaptors/CyclicFrame.h"
#include "molassembler/Temple/Adaptors/SequentialPairs.h"
#include "molassembler/Temple/constexpr/Math.h"
#include "molassembler/Temple/constexpr/Numeric.h"
#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/OrderedPair.h"
#include "molassembler/Temple/Random.h"
#include "molassembler/Temple/Stl17.h"

#include "molassembler/AngstromPositions.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/Detail/Cartesian.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/ExplicitBoundsGraph.h"
#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/Molecule.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/Stereopermutators/AbstractPermutations.h"
#include "molassembler/Stereopermutators/CycleFeasibility.h"
#include "molassembler/Stereopermutators/FeasiblePermutations.h"
#include "molassembler/Stereopermutators/ShapeVertexMaps.h"

#include <Eigen/Geometry>
#include <algorithm>

namespace Scine {
namespace Molassembler {
namespace {

template<typename ... Inds>
inline auto orderMappedSequence(
  const std::unordered_map<AtomIndex, AtomIndex>& indexMap,
  Inds ... inds
) {
  std::array<AtomIndex, sizeof...(inds)> indices {{
    indexMap.at(
      static_cast<AtomIndex>(inds)
    )...
  }};

  if(indices.front() > indices.back()) {
    std::reverse(indices.begin(), indices.end());
  }

  return indices;
}

bool piPeriodicFPCompare(const double a, const double b) {
  /* Reduces everything to [-pi, pi) bounds and then compares
   */
  auto reduceToBounds = [](const double x) -> double {
    return x - std::floor((x + M_PI)/(2 * M_PI)) * 2 * M_PI;
  };

  return std::fabs(reduceToBounds(a) - reduceToBounds(b)) < 1e-10;
}

/*
 * @note If we use tuple_element_t for SFINAE here, this works only for pair and
 * tuple. This way, it works for all types that have a first and second member,
 * of which the relevant types are std::pair and Temple::OrderedPair here.
 */
template<typename HomogeneousPair>
auto& select(
  HomogeneousPair& pair,
  bool selectFirst,
  std::enable_if_t<
    std::is_same<
      decltype(std::declval<HomogeneousPair>().first),
      decltype(std::declval<HomogeneousPair>().second)
    >::value,
    int
  >* /* enableIfPtr */ = nullptr
) {
  if(selectFirst) {
    return pair.first;
  }

  return pair.second;
}

template<typename T1, typename T2, typename ... TupleArgs>
const auto& select(
  const std::tuple<T1, T2, TupleArgs...>& tuple,
  bool selectFirst,
  std::enable_if_t<
    std::is_same<T1, T2>::value,
    int
  >* /* enableIfPtr */ = nullptr
) {
  if(selectFirst) {
    return std::get<0>(tuple);
  }

  return std::get<1>(tuple);
}

} // namespace

struct SymmetryMapHelper {
  static unsigned getSymmetryPositionOf(unsigned siteIndex, const std::vector<unsigned>& map) {
    return map.at(siteIndex);
  }

  static unsigned getSiteIndexAt(unsigned symmetryPosition, const std::vector<unsigned>& map) {
    auto findIter = std::find(
      std::begin(map),
      std::end(map),
      symmetryPosition
    );

    assert(findIter != std::end(map));

    return findIter - std::begin(map);
  }
};

std::vector<char> BondStereopermutator::Impl::charifyRankedSites_(
  const RankingInformation::RankedSitesType& sitesRanking,
  const AtomStereopermutator::ShapeMap& shapeVertexMap
) {
  std::vector<char> characters (shapeVertexMap.size());

  char currentChar = 'A';
  for(const auto& equalPrioritySet : sitesRanking) {
    for(const auto& siteIndex : equalPrioritySet) {
      characters.at(
        shapeVertexMap.at(siteIndex)
      ) = currentChar;
    }

    ++currentChar;
  }

  return characters;
}

const Stereopermutations::Composite& BondStereopermutator::Impl::composite() const {
  return composite_;
}

double BondStereopermutator::Impl::dihedral(
  const AtomStereopermutator& stereopermutatorA,
  const SiteIndex siteIndexA,
  const AtomStereopermutator& stereopermutatorB,
  const SiteIndex siteIndexB
) const {
  if(!assignment_) {
    throw std::logic_error("This stereopermutator is unassigned! Dihedrals between sites are unspecified.");
  }

  bool swapped = false;

  // Use reference wrappers to be able to swap
  using PermutatorReference = std::reference_wrapper<const AtomStereopermutator>;
  using ReferencePair = std::pair<PermutatorReference, PermutatorReference>;

  ReferencePair references {stereopermutatorA, stereopermutatorB};
  std::pair<unsigned, unsigned> siteIndices {siteIndexA, siteIndexB};
  if(stereopermutatorA.placement() == composite_.orientations().second.identifier) {
    std::swap(references.first, references.second);
    std::swap(siteIndices.first, siteIndices.second);
    swapped = true;
  }

  // Derived from possibly swapped references
  auto shapeVertexMaps = Temple::map_stl(
    references,
    [](const PermutatorReference& ref) {
      return ref.get().getShapePositionMap();
    }
  );

  std::pair<Shapes::Vertex, Shapes::Vertex> vertices;
  double dihedralAngle;

  for(const auto& dihedralTuple : composite_.dihedrals(*assignment_)) {
    std::tie(vertices.first, vertices.second, dihedralAngle) = dihedralTuple;

    if(
      shapeVertexMaps.first.indexOf(vertices.first) == siteIndices.first
      && shapeVertexMaps.second.indexOf(vertices.second) == siteIndices.second
    ) {
      return (swapped) ? -dihedralAngle : dihedralAngle;
    }
  }

  throw std::logic_error("Could not find a dihedral angle for the specified sites");
}

Stereopermutations::Composite BondStereopermutator::Impl::constructComposite_(
  const StereopermutatorList& stereopermutators,
  const BondIndex edge,
  const Alignment alignment
) {
  assert(stereopermutators.option(edge.first) && stereopermutators.option(edge.second));

  const auto& stereopermutatorA = stereopermutators.option(edge.first).value();
  const auto& stereopermutatorB = stereopermutators.option(edge.second).value();

  return {
    makeOrientationState_(stereopermutatorA, stereopermutatorB),
    makeOrientationState_(stereopermutatorB, stereopermutatorA),
    static_cast<Stereopermutations::Composite::Alignment>(alignment)
  };
}

Stereopermutations::Composite::OrientationState
BondStereopermutator::Impl::makeOrientationState_(
  const AtomStereopermutator& focalStereopermutator,
  const AtomStereopermutator& attachedStereopermutator
) {
  return {
    focalStereopermutator.getShape(),
    focalStereopermutator.getShapePositionMap().at(
      focalStereopermutator.getRanking().getSiteIndexOf(
        attachedStereopermutator.placement()
      )
    ),
    charifyRankedSites_(
      focalStereopermutator.getRanking().siteRanking,
      focalStereopermutator.getShapePositionMap()
    ),
    focalStereopermutator.placement()
  };
}

bool BondStereopermutator::Impl::cycleObviouslyInfeasible(
  const PrivateGraph& graph,
  const StereopermutatorList& stereopermutators,
  const AtomStereopermutator& firstStereopermutator,
  const AtomStereopermutator& secondStereopermutator,
  std::tuple<AtomIndex, AtomIndex, double> dihedral,
  const RankingInformation::Link& link
) {
  /* To decide if a stereopermutation is possible or not:
   * - Determine the plane into which the cycle would expand for the given
   *   dihedral angle and project the atoms of the constituting atom stereopermutators
   *   onto that plane
   * - Figure out whether the distance of the cycle atoms on a cyclic polygon
   *   to the projected atoms plus the out-of-plane distance is less than a
   *   single bond
   */
  const AtomIndex i = std::get<0>(dihedral);
  const AtomIndex j = firstStereopermutator.placement();
  const AtomIndex k = secondStereopermutator.placement();
  const AtomIndex l = std::get<1>(dihedral);
  const double phi = std::get<2>(dihedral);

  if(link.cycleSequence.size() == 3) {
    /* For triangles, modeling the cycle makes zero sense. As long as the
     * dihedral is approximately zero, the cycle must be feasible since the
     * graph is ground truth.
     */
    return std::fabs(phi) < Temple::Math::toRadians(5.0);
  }

  if(link.cycleSequence.size() == 4) {
    /* It doesn't make much more sense to model the cycle for quadrangles, since
     * AD is actually bonded and extending that distance too much is essentially
     * breaking the bond. Make a little more allowance for the quadrangle
     * than the triangle, though.
     */
    return std::fabs(phi) < Temple::Math::toRadians(20.0);
  }

  auto modelDistance = [&graph](const AtomIndex a, const AtomIndex b) -> double {
    return Bond::calculateBondDistance(
      graph.elementType(a),
      graph.elementType(b),
      graph.bondType(
        graph.edge(a, b)
      )
    );
  };

  const double a = modelDistance(i, j);
  const double b = modelDistance(j, k);
  const double c = modelDistance(k, l);

  const double alpha = DistanceGeometry::SpatialModel::siteCentralAngle(
    firstStereopermutator,
    {
      firstStereopermutator.getRanking().getSiteIndexOf(i),
      firstStereopermutator.getRanking().getSiteIndexOf(k)
    },
    graph
  );

  const double beta = DistanceGeometry::SpatialModel::siteCentralAngle(
    secondStereopermutator,
    {
      secondStereopermutator.getRanking().getSiteIndexOf(j),
      secondStereopermutator.getRanking().getSiteIndexOf(l)
    },
    graph
  );

  // Model everything in three dimensions
  const Eigen::Vector3d A = Eigen::AngleAxisd(
    alpha,
    Eigen::Vector3d::UnitZ()
  ) * Eigen::Vector3d {a, 0.0, 0.0};

  const Eigen::Vector3d B = Eigen::Vector3d::Zero();

  const Eigen::Vector3d C {b, 0.0, 0.0};

  const Eigen::Vector3d D = Eigen::AngleAxisd(
    phi,
    Eigen::Vector3d::UnitX()
  ) * (
    C + Eigen::AngleAxisd(
      M_PI - beta,
      Eigen::Vector3d::UnitZ()
    ) * Eigen::Vector3d {c, 0.0, 0.0}
  );

  const Eigen::Vector3d DMinusA = D - A;

  // P is on BC, Q on AD
  Eigen::Vector3d P, Q;

  /* If DMinusA and C are approximately parallel, getting the shortest path
   * between them is numerically dangerous
   */
  if(DMinusA.cross(C).norm() <= 1e-10) {
    // Yield P and Q from midpoints
    P = 0.5 * C;
    Q = A + 0.5 * DMinusA;
  } else {
    // Find shortest path between BC and AD
    const double mu = -(
      A.y() * DMinusA.y() + A.z() * DMinusA.z()
    ) / (
      DMinusA.y() * DMinusA.y() + DMinusA.z() * DMinusA.z()
    );

    const double lambda = (A.x() + mu * DMinusA.x()) / b;

    // Limit mu and lambda onto their line segments
    double mu_m = mu;
    if(mu <= 0) {
      mu_m = - A.x() / DMinusA.x();
    } else if(mu >= 1) {
      mu_m = (b - A.x()) / DMinusA.x();
    }

    const double lambda_m = Temple::stl17::clamp(lambda, 0.0, 1.0);

    P = lambda_m * C;
    Q = A + mu_m * DMinusA;
  }

  /* The plane the cycle expands into is now f = Plane(A, D, P). So, it's time
   * to project B and C onto the plane and then determine in-plane distances to
   * A and D.
   *
   * Plane defining vectors:
   *  m_1 = D - A = DMinusA
   *  m_2 = Q - P
   */
  const Eigen::Vector3d planeNormal = [&](){
    /* If the dihedral is an odd multiple of pi, then A, D and P become
     * collinear. The auxiliary vector used to define the plane normal is then
     * the unit z vector
     */
    if(std::fabs(std::fabs(phi) - M_PI) <= 1e-10) {
      return DMinusA.cross(Eigen::Vector3d::UnitZ()).normalized();
    }

    return DMinusA.cross(Q - P).normalized();
  }();

  // NOTE: These are signed distances!
  const double d_Bf = planeNormal.dot(A-B);
  const double d_Cf = planeNormal.dot(A-C);

  const Eigen::Vector3d projectedB = d_Bf * planeNormal + B;
  const Eigen::Vector3d projectedC = d_Cf * planeNormal + C;

  /* Create input types for cycle feasibility test */
  std::vector<Stereopermutators::BaseAtom> bases {
    // Base of B
    Stereopermutators::BaseAtom {
      graph.elementType(j),
      (A - projectedB).norm(),
      (D - projectedB).norm(),
      std::fabs(d_Bf)
    },
    // Base of C
    Stereopermutators::BaseAtom {
      graph.elementType(k),
      (A - projectedC).norm(),
      (D - projectedC).norm(),
      std::fabs(d_Cf)
    }
  };

  /* Mangle cycle sequence to get {i ... l} */
  auto cycleIndices = link.cycleSequence;
  Temple::inplace::remove_if(
    cycleIndices,
    [&](const AtomIndex x) -> bool {
      return x == j || x == k;
    }
  );
  const auto iIter = std::find(
    std::begin(cycleIndices),
    std::end(cycleIndices),
    i
  );
  assert(iIter != std::end(cycleIndices));
  std::rotate(
    std::begin(cycleIndices),
    iIter,
    std::end(cycleIndices)
  );
  if(cycleIndices.back() != l) {
    assert(cycleIndices.at(1) == l);
    std::reverse(
      std::begin(cycleIndices) + 1,
      std::end(cycleIndices)
    );
  }
  assert(cycleIndices.back() == l);

  const auto elementTypes = Temple::map(
    cycleIndices,
    [&](const AtomIndex x) -> Utils::ElementType {
      return graph.elementType(x);
    }
  );

  auto cycleEdgeLengths = Temple::map(
    Temple::adaptors::sequentialPairs(cycleIndices),
    modelDistance
  );
  cycleEdgeLengths.push_back((D-A).norm());

  // std::cout << "Cycle " << Temple::stringify(link.cycleSequence)
  //   << "(a, b, c) = (" << a << ", " << b << ", " << c << "), "
  //   << "(alpha, beta, phi) = (" << Temple::Math::toDegrees(alpha) << ", " << Temple::Math::toDegrees(beta) << ", " << Temple::Math::toDegrees(phi) << ") contradicts graph: " << fail << "\n";

  return Stereopermutators::cycleModelContradictsGraph(
    elementTypes,
    cycleEdgeLengths,
    bases
  );
}

// bool BondStereopermutator::Impl::cycleObviouslyInfeasible(
//   const PrivateGraph& graph,
//   const StereopermutatorList& stereopermutators,
//   const AtomStereopermutator& firstStereopermutator,
//   const AtomStereopermutator& secondStereopermutator,
//   std::tuple<AtomIndex, AtomIndex, double> dihedral,
//   const RankingInformation::Link& link
// ) {
//   /* Now to decide if a stereopermutation is possible or not:
//    * - Build a spatial model of each cycle including bond distances and angles
//    *   using SpatialModel's methods
//    * - Extract a BoundsList from it
//    * - Populate an ExplicitBoundsGraph
//    * - Smooth it and check for triangle inequality violations
//    */
//
//   const unsigned C = link.cycleSequence.size();
//
//   // Build a map reducing all cycle-involved atom indices to 1..C
//   std::unordered_map<AtomIndex, AtomIndex> indexReductionMap;
//   for(unsigned i = 0; i < C; ++i) {
//     indexReductionMap.emplace(
//       link.cycleSequence.at(i),
//       i
//     );
//   }
//
//   /* Build an inner graph representing only the cycle */
//   PrivateGraph minimalInner(C);
//   // Copy element types
//   for(const AtomIndex i : link.cycleSequence) {
//     minimalInner.elementType(
//       indexReductionMap.at(i)
//     ) = graph.elementType(i);
//   }
//
//   // Copy bonds
//   Temple::forEach(
//     Temple::adaptors::cyclicFrame<2>(link.cycleSequence),
//     [&](const AtomIndex i, const AtomIndex j) {
//       minimalInner.addEdge(
//         indexReductionMap.at(i),
//         indexReductionMap.at(j),
//         graph.bondType(
//           graph.edge(i, j)
//         )
//       );
//     }
//   );
//
//   /* Build a spatial model */
//   // Model all bond distances
//   DistanceGeometry::SpatialModel::BoundsMapType<2> bondBounds;
//   Temple::forEach(
//     Temple::adaptors::cyclicFrame<2>(link.cycleSequence),
//     [&](const AtomIndex i, const AtomIndex j) {
//       const BondType bondType = graph.bondType(
//         graph.edge(i, j)
//       );
//
//       double bondDistance = Bond::calculateBondDistance(
//         graph.elementType(i),
//         graph.elementType(j),
//         bondType
//       );
//
//       double absoluteVariance = bondDistance * DistanceGeometry::SpatialModel::bondRelativeVariance;
//       bondBounds.emplace(
//         orderMappedSequence(indexReductionMap, i, j),
//         DistanceGeometry::SpatialModel::makeBoundsFromCentralValue(
//           bondDistance,
//           absoluteVariance
//         )
//       );
//     }
//   );
//
//   // Model all angles
//   DistanceGeometry::SpatialModel::BoundsMapType<3> angleBounds;
//   Temple::forEach(
//     Temple::adaptors::cyclicFrame<3>(link.cycleSequence),
//     [&](const AtomIndex i, const AtomIndex j, const AtomIndex k) {
//       // Is there an AtomStereopermutator here?
//       DistanceGeometry::ValueBounds angleValueBounds {
//         Shapes::smallestAngle,
//         M_PI
//       };
//
//       if(auto permutatorOption = stereopermutators.option(j)) {
//         const RankingInformation& ranking = permutatorOption->getRanking();
//         // Make sure this is as simple a case as possible to model
//         unsigned siteOfI = ranking.getSiteIndexOf(i);
//         unsigned siteOfK = ranking.getSiteIndexOf(k);
//
//         // Both sites are single-index, not haptic
//         if(
//           ranking.sites.at(siteOfI).size() == 1
//           && ranking.sites.at(siteOfK).size() == 1
//         ) {
//           if(permutatorOption->assigned()) {
//             angleValueBounds = DistanceGeometry::SpatialModel::makeBoundsFromCentralValue(
//               permutatorOption->angle(siteOfI, siteOfK),
//               DistanceGeometry::SpatialModel::angleAbsoluteVariance
//             );
//           } else {
//             angleValueBounds = {
//               Shapes::minimumAngle(permutatorOption->getShape()),
//               Shapes::maximumAngle(permutatorOption->getShape())
//             };
//           }
//         }
//       }
//
//       angleBounds.emplace(
//         orderMappedSequence(indexReductionMap, i, j, k),
//         angleValueBounds
//       );
//     }
//   );
//
//   // Model dihedrals
//   DistanceGeometry::SpatialModel::BoundsMapType<4> dihedralBounds;
//   const double dihedralAngle = std::get<2>(dihedral);
//   auto dihedralValueBounds = DistanceGeometry::SpatialModel::makeBoundsFromCentralValue(
//     dihedralAngle,
//     DistanceGeometry::SpatialModel::dihedralAbsoluteVariance
//   );
//
//   const AtomIndex i = std::get<0>(dihedral);
//   const AtomIndex j = firstStereopermutator.placement();
//   const AtomIndex k = secondStereopermutator.placement();
//   const AtomIndex l = std::get<1>(dihedral);
//
//   assert(Temple::makeContainsPredicate(link.cycleSequence)(i));
//   assert(Temple::makeContainsPredicate(link.cycleSequence)(l));
//
//   dihedralBounds.emplace(
//     orderMappedSequence(indexReductionMap, i, j, k, l),
//     dihedralValueBounds
//   );
//
//   // Create pairwise bounds from internal coordinates
//   auto pairwiseBounds = DistanceGeometry::SpatialModel::makePairwiseBounds(
//     C,
//     DistanceGeometry::SpatialModel::BoundsMapType<2> {},
//     bondBounds,
//     angleBounds,
//     dihedralBounds
//   );
//
//   // Model in a bounds matrix
//   DistanceGeometry::ExplicitBoundsGraph boundsGraph {
//     minimalInner,
//     pairwiseBounds
//   };
//
//   // Try to smooth the distance bounds with triangle inequalities
//   auto distanceBoundsResult = boundsGraph.makeDistanceBounds();
//
//   /*if(distanceBoundsResult) {
//     std::cout << "Link " << Temple::stringify(link.cycleSequence) << " is viable on " << j << ", " << k << " with dihedral " << dihedralAngle << "\n";
//     std::cout << pairwiseBounds << "\n";
//     std::cout << distanceBoundsResult.value() << "\n";
//
//     auto distanceMatrix = boundsGraph.makeDistanceMatrix(randomnessEngine());
//     if(distanceMatrix) {
//       std::cout << distanceMatrix.value() << "\n";
//     }
//
//     std::cout << "\n";
//   }*/
//
//   return !(distanceBoundsResult.operator bool());
// }

std::vector<unsigned> BondStereopermutator::Impl::notObviouslyInfeasibleStereopermutations(
  const PrivateGraph& graph,
  const StereopermutatorList& stereopermutators,
  const Stereopermutations::Composite& composite
) {
  const unsigned compositePermutations = composite.permutations();

  auto permutatorReferences = composite.orientations().map(
    [&](const auto& orientationState) -> const AtomStereopermutator& {
      auto permutatorOption = stereopermutators.option(orientationState.identifier);
      if(!permutatorOption) {
        throw std::logic_error("Atom stereopermutators are not present for this composite!");
      }
      return *permutatorOption;
    }
  );

  /* Before iterating through all possible stereopermutations, figure out
   * if substituents of both stereopermutators are fused somehow, i.e. if there
   * are cycles involving the bond between A and B.
   */
  auto links = GraphAlgorithms::siteLinks(
    graph,
    permutatorReferences.first,
    permutatorReferences.second
  );

  /* If there are no links between the substituents of the stereopermutators,
   * all permutations are presumably feasible.
   */
  if(links.empty()) {
    return Temple::iota<unsigned>(compositePermutations);
  }

  /* Try to match the shape position of the composite dihedral to the link's
   * site index pair using the atom stereopermutator shape position maps
   *
   * Note that there may be links that are not involved in dihedrals (i.e.
   * they are not part of the group of shape positions that interact in the
   * dihedral, e.g. the position trans in a triangle - square combination.)
   */
  auto getDihedralInformation = [&](
    const std::vector<Stereopermutations::Composite::DihedralTuple>& dihedrals,
    const RankingInformation::Link& link
  ) -> boost::optional<std::tuple<AtomIndex, AtomIndex, double>> {
    const auto& firstSiteIndices = permutatorReferences.first.getRanking().sites.at(link.sites.first);
    const auto& secondSiteIndices = permutatorReferences.second.getRanking().sites.at(link.sites.second);

    // We can't decide a dihedral for haptic sites here
    if(firstSiteIndices.size() > 1 || secondSiteIndices.size() > 1) {
      return boost::none;
    }

    const std::pair<unsigned, unsigned> linkShapePositions {
      permutatorReferences.first.getShapePositionMap().at(link.sites.first),
      permutatorReferences.second.getShapePositionMap().at(link.sites.second)
    };

    // Look for a composite dihedral matching the shape positions
    auto findIter = std::find_if(
      std::begin(dihedrals),
      std::end(dihedrals),
      [&](const auto& dihedral) -> bool {
        return (
          std::get<0>(dihedral) == std::get<0>(linkShapePositions)
          && std::get<1>(dihedral) == std::get<1>(linkShapePositions)
        );
      }
    );

    // There is no matching dihedral for this link (this is okay!)
    if(findIter == std::end(dihedrals)) {
      return boost::none;
    }

    const double dihedralAngle = std::get<2>(*findIter);
    return std::make_tuple(
      firstSiteIndices.front(),
      secondSiteIndices.front(),
      dihedralAngle
    );
  };

  std::vector<unsigned> viableStereopermutations;
  for(
    unsigned stereopermutationIndex = 0;
    stereopermutationIndex < compositePermutations;
    ++stereopermutationIndex
  ) {
    if(
      !Temple::any_of(
        links,
        [&](const RankingInformation::Link& link) -> bool {
          auto dihedralInformationOption = getDihedralInformation(
            composite.dihedrals(stereopermutationIndex),
            link
          );

          if(!dihedralInformationOption) {
            /* Haptic sites involved, we don't deal with those currently, we
             * blanket accept them. Or the link is irrelevant to the dihedral,
             * so it doesn't matter to the viability of the stereopermutation.
             */
            return false;
          }

          return cycleObviouslyInfeasible(
            graph,
            stereopermutators,
            permutatorReferences.first,
            permutatorReferences.second,
            *dihedralInformationOption,
            link
          );
        }
      )
    ) {
      viableStereopermutations.push_back(stereopermutationIndex);
    }
  }

  return viableStereopermutations;
}

/* Constructors */
BondStereopermutator::Impl::Impl(
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const BondIndex edge,
  Alignment alignment
) : composite_ {
      makeOrientationState_(stereopermutatorA, stereopermutatorB),
      makeOrientationState_(stereopermutatorB, stereopermutatorA),
      static_cast<Stereopermutations::Composite::Alignment>(alignment)
    },
    edge_(edge),
    feasiblePermutations_(
      Temple::iota<unsigned>(composite_.permutations())
    ),
    assignment_(boost::none)
{}

BondStereopermutator::Impl::Impl(
  const PrivateGraph& graph,
  const StereopermutatorList& stereopermutators,
  const BondIndex edge,
  Alignment alignment
) : composite_(
      constructComposite_(stereopermutators, edge, alignment)
    ),
    edge_(edge),
    feasiblePermutations_(
      notObviouslyInfeasibleStereopermutations(
        graph, stereopermutators, composite_
      )
    ),
    assignment_(boost::none)
{}

/* Public members */
/* Modification */
void BondStereopermutator::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment && assignment.value() >= numAssignments()) {
    /* The distinction here between numAssignments and composite_.permutations()
     * is important because if the composite is isotropic, we simulate that
     * there is only a singular assignment (see numAssignments), although
     * all permutations are generated and present for fitting anyway.
     *
     * If this were composite_.permutations(), we would accept assignments other
     * than zero if the underlying composite is isotropic, but not yield the
     * same assignment index when asked which assignment is currently set
     * in assigned().
     */
    throw std::out_of_range("Supplied assignment index is out of range");
  }

  assignment_ = assignment;
}

void BondStereopermutator::Impl::assignRandom(Random::Engine& engine) {
  const unsigned A = numAssignments();
  if(A == 0) {
    throw std::logic_error("Cannot randomly assign a stereopermutator without feasible stereopermutations");
  }

  if(A == 1) {
    assign(0);
  } else {
    assign(
      Temple::Random::getSingle<unsigned>(0, A - 1, engine)
    );
  }
}

void BondStereopermutator::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
  /* Composite's OrientationState identifiers change (these have no impact on
   * assignment ordering however, so that should be safe)
   */
  composite_.applyIdentifierPermutation(permutation);

  // Edge we sit on changes according to the permutation
  edge_ = BondIndex {
    permutation.at(edge_.first),
    permutation.at(edge_.second)
  };

  // Assignment does not change
}

void BondStereopermutator::Impl::fit(
  const AngstromPositions& angstromWrapper,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const FittingMode mode
) {
  assert(stereopermutatorA.placement() != stereopermutatorB.placement());

  // Early exit
  if(composite_.permutations() == 0) {
    assignment_ = boost::none;
    return;
  }

  // Can the following selections be done with a single branch?
  const AtomStereopermutator& firstStereopermutator = (
    stereopermutatorA.placement() == composite_.orientations().first.identifier
    ? stereopermutatorA
    : stereopermutatorB
  );

  const AtomStereopermutator& secondStereopermutator = (
    stereopermutatorB.placement() == composite_.orientations().second.identifier
    ? stereopermutatorB
    : stereopermutatorA
  );

  auto makeSitePositions = [&angstromWrapper](const AtomStereopermutator& permutator) -> Eigen::Matrix<double, 3, Eigen::Dynamic> {
    const unsigned S = permutator.getRanking().sites.size();
    assert(S == Shapes::size(permutator.getShape()));
    Eigen::Matrix<double, 3, Eigen::Dynamic> sitePositions(3, S);
    for(unsigned i = 0; i < S; ++i) {
      sitePositions.col(i) = cartesian::averagePosition(
        angstromWrapper.positions,
        permutator.getRanking().sites.at(i)
      );
    }
    return sitePositions;
  };

  // For all atoms making up a site, decide on the spatial average position
  auto firstSitePositions = makeSitePositions(firstStereopermutator);
  auto secondSitePositions = makeSitePositions(secondStereopermutator);

  Shapes::Vertex firstShapeVertex;
  Shapes::Vertex secondShapeVertex;
  double dihedralAngle;

  double bestPenalty = std::numeric_limits<double>::max();
  std::vector<unsigned> bestStereopermutation;

  for(unsigned feasiblePermutationIndex : feasiblePermutations_) {
    double penalty = 0.0;
    for(const auto& dihedralTuple : composite_.dihedrals(feasiblePermutationIndex)) {
      std::tie(firstShapeVertex, secondShapeVertex, dihedralAngle) = dihedralTuple;

      // Get site index of leftSymmetryPosition in left
      const SiteIndex firstSite = firstStereopermutator.getShapePositionMap().indexOf(firstShapeVertex);
      const SiteIndex secondSite = secondStereopermutator.getShapePositionMap().indexOf(secondShapeVertex);

      /* Dihedral angle differences aren't as easy as |b - a|, since
       * dihedrals are defined over (-pi, pi], so in the worst case
       *
       *   a = pi, b = -pi + delta
       *   |(-pi + delta) - pi| = 2 pi - delta.
       *
       * Any difference has to be modified by a 2 pi-periodic function
       * calculating the distance to the nearest multiple of pi (like a
       * triangle wave function) or we have to ensure that the difference is
       * also in the definition interval and then calculate the distance to
       * zero.
       */

      const double measuredDihedral = cartesian::dihedral(
        firstSitePositions.col(firstSite),
        angstromWrapper.positions.row(firstStereopermutator.placement()),
        angstromWrapper.positions.row(secondStereopermutator.placement()),
        secondSitePositions.col(secondSite)
      );

      double dihedralDifference = measuredDihedral - dihedralAngle;

      // + pi is part of the definition interval, so use greater than
      if(dihedralDifference > M_PI) {
        dihedralDifference -= 2 * M_PI;
      }

      // - pi is not part of the definition interval, so use less than or equal
      if(dihedralDifference <= -M_PI) {
        dihedralDifference += 2 * M_PI;
      }

      penalty += std::fabs(dihedralDifference);
    }

    /* Check if this penalty is within acceptance threshold if fitting mode
     * encompasses this
     *
     * The logic here is that the acceptable deviation per dihedral should
     * depend on the order of the Composite (i.e. if there are three
     * substituents at one side, the deviation per dihedral should be smaller
     * than if there were only two).
     */
    if(
      mode == FittingMode::Thresholded
      && (
        penalty / composite_.dihedrals(feasiblePermutationIndex).size()
        > assignmentAcceptanceParameter * 2 * M_PI / composite_.order()
      )
    ) {
      continue;
    }

    if(penalty < bestPenalty) {
      bestPenalty = penalty;
      bestStereopermutation = {feasiblePermutationIndex};
    } else if(penalty == bestPenalty) {
      bestStereopermutation.push_back(feasiblePermutationIndex);
    }
  }

  /* The best stereopermutation must be singular, no other may match it in
   * fit value, otherwise assignment is ambiguous
   */
  if(bestStereopermutation.size() == 1) {
    /* Retrieve the index of the best stereopermutation from among the feasible
     * permutations
     */
    const unsigned stereopermutation = bestStereopermutation.front();
    assignment_ = (
      Temple::find(feasiblePermutations_, stereopermutation)
      - std::begin(feasiblePermutations_)
    );
  } else {
    // On ambiguous matching, dis-assign the stereopermutator
    assignment_ = boost::none;
  }
}

void BondStereopermutator::Impl::propagateGraphChange(
  const AtomStereopermutatorPropagatedState& oldPermutatorState,
  const AtomStereopermutator& newPermutator,
  const PrivateGraph& graph,
  const StereopermutatorList& permutators
) {
  const RankingInformation& oldRanking = std::get<0>(oldPermutatorState);
  const Stereopermutators::Abstract& oldAbstract = std::get<1>(oldPermutatorState);
  const Stereopermutators::Feasible& oldFeasible = std::get<2>(oldPermutatorState);
  const boost::optional<unsigned>& oldAssignment = std::get<3>(oldPermutatorState);

  // We assume that the supplied permutators (or their state) were assigned
  assert(oldAssignment);
  assert(newPermutator.assigned());

  /* We assume the old and new symmetry are of the same size (i.e. this is
   * a ranking change propagation, not a substituent addition / removal
   * propagation)
   */
  assert(oldRanking.sites.size() == newPermutator.getRanking().sites.size());

  using OrientationState = Stereopermutations::Composite::OrientationState;

  /* Construct a new Composite with the new information */
  bool changedIsFirstInOldOrientations = (
    composite_.orientations().first.identifier == newPermutator.placement()
  );

  const OrientationState& oldOrientation = select(
    composite_.orientations(),
    changedIsFirstInOldOrientations
  );

  // Reuse the OrientationState of the "other" atom stereopermutator
  const OrientationState& unchangedOrientation = select(
    composite_.orientations(),
    !changedIsFirstInOldOrientations
  );

  // Generate a new OrientationState for the modified stereopermutator
  Stereopermutations::Composite::OrientationState possiblyModifiedOrientation {
    newPermutator.getShape(),
    newPermutator.getShapePositionMap().at(
      newPermutator.getRanking().getSiteIndexOf(
        unchangedOrientation.identifier
      )
    ),
    charifyRankedSites_(
      newPermutator.getRanking().siteRanking,
      newPermutator.getShapePositionMap()
    ),
    newPermutator.placement()
  };

  // In case nothing has changed, we are done and can stop
  if(oldOrientation == possiblyModifiedOrientation) {
    return;
  }

  /* Generate a new set of permutations (note this may reorder its
   * arguments into the orientations() member)
   */
  Stereopermutations::Composite newComposite {
    possiblyModifiedOrientation,
    unchangedOrientation
  };

  // feasibility has to be rechecked
  std::vector<unsigned> newFeasiblePermutations;
  newFeasiblePermutations = notObviouslyInfeasibleStereopermutations(
    graph,
    permutators,
    newComposite
  );

  /* If the new composite is isotropic, there is no reason to try to find
   * an assignment, we can choose any.
   */
  if(newComposite.isIsotropic()) {
    composite_ = std::move(newComposite);
    feasiblePermutations_ = std::move(newFeasiblePermutations);
    if(!feasiblePermutations_.empty()) {
      assignment_ = 0;
    }
    return;
  }

  /* If this BondStereopermutator is unassigned, and the permutator is not
   * isotropic after the ranking change, then this permutator stays unassigned.
   */
  if(assignment_ == boost::none) {
    composite_ = std::move(newComposite);
    feasiblePermutations_ = std::move(newFeasiblePermutations);
    return;
  }

  /* Find the old permutation in the set of new permutations
   * - Since composite only offers dihedral information in terms of symmetry
   *   positions, we have to translate these back into site indices
   * - Transform symmetry positions through the sequence
   *
   *   old symmetry position
   *   -> site index (using old permutation state)
   *   -> atom index (using old ranking)
   *   -> site index (using new ranking)
   *   -> new symmetry position (using new permutation state)
   *
   *   and then compare dihedral values between the composites. This scheme
   *   has works even if the fused position has changed.
   * - Since newComposite may have reordered the OrientationStates, we have
   *   to be careful which part of the DihedralTuple's we extract symmetry
   *   positions from.
   * - Inversions of the dihedral defining sequence do not invert the sign of
   *   the dihedral
   */
  bool modifiedOrientationIsFirstInNewComposite = (
    newComposite.orientations().first.identifier
    == possiblyModifiedOrientation.identifier
  );

  using DihedralTuple = Stereopermutations::Composite::DihedralTuple;
  // This permutator is assigned since that is ensured a few lines earlier
  const std::vector<DihedralTuple>& oldDihedralList = composite_.dihedrals(
    assignment_.value()
  );

  auto oldSymmetryPositionToSiteMap = shapeVertexToSiteIndexMap(
    oldAbstract.permutations.list.at(
      oldFeasible.indices.at(
        oldAssignment.value()
      )
    ),
    oldAbstract.canonicalSites,
    oldRanking.links
  );

  auto getNewShapeVertex = [&](Shapes::Vertex oldVertex) -> Shapes::Vertex {
    const SiteIndex oldSiteIndex = oldSymmetryPositionToSiteMap.at(oldVertex);

    const std::vector<AtomIndex>& oldSite = oldRanking.sites.at(oldSiteIndex);

    // We assume here that atom indices making up sites are sorted
    assert(std::is_sorted(std::begin(oldSite), std::end(oldSite)));

    /* Find this site in the new ranking
     * - In case there is truly only a ranking change (not a rearrangement in
     *   terms of haptic sites or so), then there should be a vector with
     *   exactly the same elements.
     */
    const auto& newRankingSites = newPermutator.getRanking().sites;
    const auto findSiteIter = std::find(
      std::begin(newRankingSites),
      std::end(newRankingSites),
      oldSite
    );

    assert(findSiteIter != std::end(newRankingSites));

    const SiteIndex newSiteIndex = SiteIndex(findSiteIter - std::begin(newRankingSites));
    return newPermutator.getShapePositionMap().at(newSiteIndex);
  };

  // Map the set of old dihedral tuples into the new space
  std::vector<DihedralTuple> newCompositeDihedrals;
  newCompositeDihedrals.reserve(oldDihedralList.size());
  for(const DihedralTuple& oldDihedral : oldDihedralList) {
    const Shapes::Vertex changedVertex = select(
      oldDihedral,
      changedIsFirstInOldOrientations
    );

    const Shapes::Vertex unchangedVertex = select(
      oldDihedral,
      !changedIsFirstInOldOrientations
    );

    const Shapes::Vertex newShapeVertex = getNewShapeVertex(changedVertex);

    if(modifiedOrientationIsFirstInNewComposite) {
      newCompositeDihedrals.emplace_back(
        newShapeVertex,
        unchangedVertex,
        std::get<2>(oldDihedral)
      );
    } else {
      newCompositeDihedrals.emplace_back(
        unchangedVertex,
        newShapeVertex,
        std::get<2>(oldDihedral)
      );
    }
  }

  /* Sort the dihedrals for easy comparison. Composite already generates its
   * dihedrals in a sorted fashion, so newComposite's individual permutations'
   * dihedrals need not be sorted.
   */
  std::sort(
    std::begin(newCompositeDihedrals),
    std::end(newCompositeDihedrals)
  );

  /* Find the matching stereopermutation by finding a bijective mapping from
   * all old dihedral tuples of the old stereopermutation to all dihedral
   * tuples of a stereopermutation within the new Composite. Since all lists of
   * DihedralTuples are sorted, we can use the lexicographical vector
   * comparison operator that in turn calls the lexicographical tuple
   * comparison operator:
   *
   * However, the match may not be floating-point exact, and can have pi
   * periodicities!
   */
  auto matchIter = std::find_if(
    std::begin(newComposite),
    std::end(newComposite),
    [&](const std::vector<DihedralTuple>& dihedrals) -> bool {
      return std::lexicographical_compare(
        std::begin(newCompositeDihedrals),
        std::end(newCompositeDihedrals),
        std::begin(dihedrals),
        std::end(dihedrals),
        [](const DihedralTuple& a, const DihedralTuple& b) -> bool {
          return (
            std::get<0>(a) == std::get<0>(b)
            && std::get<1>(a) == std::get<1>(b)
            && piPeriodicFPCompare(std::get<2>(a), std::get<2>(b))
          );
        }
      );
    }
  );

  // There should always be a match
  if(matchIter == std::end(newComposite)) {
    throw std::logic_error("Bug: no match found in new composite.");
  }

  // Overwrite class state
  assignment_ = matchIter - std::begin(newComposite);
  composite_ = std::move(newComposite);
  feasiblePermutations_ = std::move(newFeasiblePermutations);
}

/* Information */
BondStereopermutator::Alignment BondStereopermutator::Impl::alignment() const {
  return static_cast<BondStereopermutator::Alignment>(
    composite_.alignment()
  );
}

boost::optional<unsigned> BondStereopermutator::Impl::assigned() const {
  /* If the underlying composite is isotropic, it does not matter which of those
   * permutations by shape position is the factual spatial arrangement (since
   * they are all rotationally equivalent). We have to spoof that there is only
   * one arrangement in this case (although we need all of them for spatial
   * fitting).
   */
  if(assignment_ && (composite_.isIsotropic() && !feasiblePermutations_.empty())) {
    return 0U;
  }

  return assignment_;
}

bool BondStereopermutator::Impl::hasSameCompositeOrientation(const BondStereopermutator::Impl& other) const {
  return composite_ == other.composite_;
}

boost::optional<unsigned> BondStereopermutator::Impl::indexOfPermutation() const {
  if(assignment_ && composite_.isIsotropic()) {
    return 0U;
  }

  if(!assignment_) {
    return boost::none;
  }

  return feasiblePermutations_.at(*assignment_);
}

unsigned BondStereopermutator::Impl::numAssignments() const {
  if(composite_.isIsotropic() && !feasiblePermutations_.empty()) {
    return 1;
  }

  return feasiblePermutations_.size();
}

unsigned BondStereopermutator::Impl::numStereopermutations() const {
  if(composite_.isIsotropic()) {
    return 1;
  }

  return composite_.permutations();
}

std::string BondStereopermutator::Impl::info() const {
  using namespace std::string_literals;

  std::string returnString =  "";

  returnString += std::to_string(composite_.orientations().first.identifier);
  returnString += "-";
  returnString += std::to_string(composite_.orientations().second.identifier);

  const unsigned A = numAssignments();

  if(A == 1) {
    returnString += ": Is non-stereogenic.";
  } else {
    returnString += ": Is ";
    if(assignment_) {
      returnString += std::to_string(assignment_.value());
    } else {
      returnString += "u";
    }

    returnString += " ("s + std::to_string(A);

    const unsigned P = numStereopermutations();
    if(P != A) {
      returnString += ", "s + std::to_string(P);
    }
    returnString += ")";
  }

  return returnString;
}

std::string BondStereopermutator::Impl::rankInfo() const {
  using namespace std::string_literals;

  return (
    "B-"s + std::to_string(numStereopermutations())
    + "-"s + (
      assigned()
      ? std::to_string(assigned().value())
      : "u"s
    )
  );
}

BondIndex BondStereopermutator::Impl::placement() const {
  // Return a standard form of smaller first
  return edge_;
}

} // namespace Molassembler
} // namespace Scine
