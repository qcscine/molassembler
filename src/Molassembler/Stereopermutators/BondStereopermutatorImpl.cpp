/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutators/BondStereopermutatorImpl.h"

#include "Molassembler/Shapes/Data.h"

#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Adaptors/CyclicFrame.h"
#include "Molassembler/Temple/Adaptors/SequentialPairs.h"
#include "Molassembler/Temple/constexpr/Math.h"
#include "Molassembler/Temple/constexpr/Numeric.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/OrderedPair.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stl17.h"
#include "Molassembler/Temple/Optionals.h"

#include "Molassembler/AngstromPositions.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/Detail/Cartesian.h"
#include "Molassembler/DistanceGeometry/SpatialModel.h"
#include "Molassembler/DistanceGeometry/ExplicitBoundsGraph.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Modeling/BondDistance.h"
#include "Molassembler/Modeling/CommonTrig.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/RankingInformation.h"
#include "Molassembler/Stereopermutators/AbstractPermutations.h"
#include "Molassembler/Stereopermutators/CycleFeasibility.h"
#include "Molassembler/Stereopermutators/FeasiblePermutations.h"
#include "Molassembler/Stereopermutators/ShapeVertexMaps.h"

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

std::pair<BondStereopermutator::FittingReferences, BondStereopermutator::FittingReferences> align(
  std::pair<BondStereopermutator::FittingReferences, BondStereopermutator::FittingReferences> references,
  const Stereopermutations::Composite& composite
) {
  const auto orientationIdentifiers = composite.orientations().map(
    [](const auto& x) { return x.identifier; }
  );

  const auto placements = Temple::map(
    references,
    [](const auto& r) -> AtomIndex { return r.stereopermutator.placement(); }
  );
  assert(placements.first != placements.second);

  if(orientationIdentifiers == placements) {
    return references;
  }

  const auto reversedPlacements = std::make_pair(placements.second, placements.first);
  if(orientationIdentifiers == reversedPlacements) {
    return std::make_pair(
      references.second,
      references.first
    );
  }

  throw std::runtime_error("Fitting references and composite orientation identifiers are mismatched!");
}

template<typename T, typename U, typename F>
auto zipMapPair(
  const T& t,
  const U& u,
  F&& f
) {
  return std::make_pair(
    f(t.first, u.first),
    f(t.second, u.second)
  );
}

} // namespace

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
  auto shapeVertexMaps = Temple::map(
    references,
    [](const PermutatorReference& ref) {
      return ref.get().getShapePositionMap();
    }
  );

  std::pair<Shapes::Vertex, Shapes::Vertex> vertices;
  double dihedralAngle;

  for(const auto& dihedralTuple : composite_.allPermutations().at(*assignment_).dihedrals) {
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
    makeOrientationState_(stereopermutatorA, stereopermutatorA.getShapePositionMap(), stereopermutatorB),
    makeOrientationState_(stereopermutatorB, stereopermutatorB.getShapePositionMap(), stereopermutatorA),
    static_cast<Stereopermutations::Composite::Alignment>(alignment)
  };
}

Stereopermutations::Composite::OrientationState
BondStereopermutator::Impl::makeOrientationState_(
  const AtomStereopermutator& focalStereopermutator,
  const AtomStereopermutator::ShapeMap& focalShapeMap,
  const AtomStereopermutator& attachedStereopermutator
) {
  return {
    focalStereopermutator.getShape(),
    focalShapeMap.at(
      focalStereopermutator.getRanking().getSiteIndexOf(
        attachedStereopermutator.placement()
      )
    ),
    charifyRankedSites_(
      focalStereopermutator.getRanking().siteRanking,
      focalShapeMap
    ),
    focalStereopermutator.placement()
  };
}

bool BondStereopermutator::Impl::cycleObviouslyInfeasible(
  const PrivateGraph& graph,
  const StereopermutatorList& /* stereopermutators */,
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
  Eigen::Vector3d P;
  Eigen::Vector3d Q;

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

    const double lambda_m = Temple::Stl17::clamp(lambda, 0.0, 1.0);

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
  Temple::remove_if(
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
    Temple::Adaptors::sequentialPairs(cycleIndices),
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
//     Temple::Adaptors::cyclicFrame<2>(link.cycleSequence),
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
//     Temple::Adaptors::cyclicFrame<2>(link.cycleSequence),
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
//     Temple::Adaptors::cyclicFrame<3>(link.cycleSequence),
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
  const auto permutatorReferences = composite.orientations().map(
    [&](const auto& orientationState) -> const AtomStereopermutator& {
      return stereopermutators.at(orientationState.identifier);
    }
  );

  /* Before iterating through all possible stereopermutations, figure out
   * if substituents of both stereopermutators are fused somehow, i.e. if there
   * are cycles involving the bond between A and B.
   */
  const auto links = GraphAlgorithms::siteLinks(
    graph,
    permutatorReferences.first,
    permutatorReferences.second
  );

  /* If there are no links between the substituents of the stereopermutators,
   * all permutations are presumably feasible.
   */
  if(links.empty()) {
    return composite.nonEquivalentPermutationIndices();
  }

  /* Try to match the shape vertex of the composite dihedral to the link's
   * site index pair using the atom stereopermutator shape position maps
   *
   * Note that there may be links that are not involved in dihedrals (i.e.
   * they are not part of the group of shape positions that interact in the
   * dihedral, e.g. the position trans in a triangle - square combination.)
   */
  auto getDihedralInformation = [&](
    const std::vector<Stereopermutations::Composite::Permutation::DihedralTuple>& dihedrals,
    const RankingInformation::Link& link
  ) -> boost::optional<std::tuple<AtomIndex, AtomIndex, double>> {
    const auto& firstSiteIndices = permutatorReferences.first.getRanking().sites.at(link.sites.first);
    const auto& secondSiteIndices = permutatorReferences.second.getRanking().sites.at(link.sites.second);

    // We can't decide a dihedral for haptic sites here
    if(firstSiteIndices.size() > 1 || secondSiteIndices.size() > 1) {
      return boost::none;
    }

    const std::pair<Shapes::Vertex, Shapes::Vertex> linkVertices {
      permutatorReferences.first.getShapePositionMap().at(link.sites.first),
      permutatorReferences.second.getShapePositionMap().at(link.sites.second)
    };

    // Look for a composite dihedral matching the shape positions
    const auto findIter = Temple::find_if(
      dihedrals,
      [&](const auto& dihedralTuple) -> bool {
        return (
          std::get<0>(dihedralTuple) == std::get<0>(linkVertices)
          && std::get<1>(dihedralTuple) == std::get<1>(linkVertices)
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

  std::set<unsigned> viableStereopermutations;
  const auto permEnd = std::end(composite);
  for(auto permIter = std::begin(composite); permIter != permEnd; ++permIter) {
    if(permIter->rankingEquivalentTo) {
      continue;
    }

    if(
      !Temple::any_of(
        links,
        [&](const RankingInformation::Link& link) -> bool {
          auto infoOption = getDihedralInformation(permIter->dihedrals, link);

          if(!infoOption) {
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
            *infoOption,
            link
          );
        }
      )
    ) {
      viableStereopermutations.insert(
        composite.rankingEquivalentBase(permIter - std::begin(composite))
      );
    }
  }

  return {std::begin(viableStereopermutations), std::end(viableStereopermutations)};
}

/* Constructors */
BondStereopermutator::Impl::Impl(
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const BondIndex edge,
  Alignment alignment
) : composite_ {
      makeOrientationState_(stereopermutatorA, stereopermutatorA.getShapePositionMap(), stereopermutatorB),
      makeOrientationState_(stereopermutatorB, stereopermutatorB.getShapePositionMap(), stereopermutatorA),
      static_cast<Stereopermutations::Composite::Alignment>(alignment)
    },
    edge_(edge),
    feasiblePermutations_(composite_.nonEquivalentPermutationIndices()),
    assignment_(boost::none)
{}

BondStereopermutator::Impl::Impl(
  const PrivateGraph& graph,
  const StereopermutatorList& stereopermutators,
  const BondIndex edge,
  Alignment alignment
) : composite_(constructComposite_(stereopermutators, edge, alignment)),
    edge_(edge),
    feasiblePermutations_(
      notObviouslyInfeasibleStereopermutations(graph, stereopermutators, composite_)
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
  const SitePositionsPair& sitePositions,
  std::pair<FittingReferences, FittingReferences> fittingReferences,
  const FittingMode mode
) {
  // Early exit
  if(composite_.countNonEquivalentPermutations() == 0) {
    assignment_ = boost::none;
    return;
  }

  const auto alignedReferences = align(fittingReferences, composite_);
  Stereopermutations::Composite matchedComposite {
    makeOrientationState_(
      alignedReferences.first.stereopermutator,
      alignedReferences.first.shapeMap,
      alignedReferences.second.stereopermutator
    ),
    makeOrientationState_(
      alignedReferences.second.stereopermutator,
      alignedReferences.second.shapeMap,
      alignedReferences.first.stereopermutator
    ),
    composite_.alignment()
  };

  std::pair<Shapes::Vertex, Shapes::Vertex> shapeVertices;
  double dihedralAngle;

  double bestMisalignment = std::numeric_limits<double>::max();
  std::vector<unsigned> bestStereopermutations;

  // Expand the feasible permutations by ranking equivalent rotations
  std::unordered_map<unsigned, unsigned> permutationMapToFeasibleBase;
  for(const unsigned i : feasiblePermutations_) {
    permutationMapToFeasibleBase.emplace(i, i);
  }
  for(unsigned i = 0; i < composite_.allPermutations().size(); ++i) {
    unsigned base = composite_.rankingEquivalentBase(i);

    if(permutationMapToFeasibleBase.count(base) > 0) {
      permutationMapToFeasibleBase.emplace(i, base);
    }
  }

  for(auto mapPair : permutationMapToFeasibleBase) {
    const unsigned feasiblePermutationIndex = mapPair.first;

    double misalignment = 0.0;
    for(
      const auto& dihedralTuple :
      matchedComposite.allPermutations().at(feasiblePermutationIndex).dihedrals
    ) {
      std::tie(shapeVertices.first, shapeVertices.second, dihedralAngle) = dihedralTuple;

      const auto siteIndices = zipMapPair(
        alignedReferences,
        shapeVertices,
        [](const FittingReferences& ref, const Shapes::Vertex v) -> SiteIndex {
          return ref.shapeMap.indexOf(v);
        }
      );

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

      const double measuredDihedral = Cartesian::dihedral(
        sitePositions.first.col(siteIndices.first),
        sitePositions.first.rightCols<1>(),
        sitePositions.second.rightCols<1>(),
        sitePositions.second.col(siteIndices.second)
      );

      misalignment += Cartesian::dihedralDifference(measuredDihedral, dihedralAngle);
    }

    /* Check if this misalignment is within acceptance threshold if fitting mode
     * is thresholded
     *
     * The logic here is that the acceptable deviation per dihedral should
     * depend on the order of the Composite (i.e. if there are three
     * substituents at one side, the deviation per dihedral should be smaller
     * than if there were only two).
     */
    const double misalignmentPerDihedral = misalignment / composite_.allPermutations().at(feasiblePermutationIndex).dihedrals.size();
    const double acceptableMisalignment = assignmentAcceptanceParameter * 2 * M_PI / composite_.order();
    if(
      mode == FittingMode::Thresholded
      && misalignmentPerDihedral > acceptableMisalignment
    ) {
      continue;
    }

    if(misalignment < bestMisalignment) {
      bestMisalignment = misalignment;
      bestStereopermutations = {feasiblePermutationIndex};
    } else if(misalignment == bestMisalignment) {
      bestStereopermutations.push_back(feasiblePermutationIndex);
    }
  }

  /* The best permutation must be singular, no other may match it in fit value,
   * otherwise assignment is ambiguous
   */
  if(bestStereopermutations.size() == 1) {
    /* Retrieve the index of the best permutation from among the feasible
     * permutations
     */
    const unsigned feasiblePermutationIndex = permutationMapToFeasibleBase.at(bestStereopermutations.front());
    auto findIter = Temple::find(feasiblePermutations_, feasiblePermutationIndex);
    assert(findIter != std::end(feasiblePermutations_));
    assignment_ = findIter - std::begin(feasiblePermutations_);
  } else {
    // On ambiguous matching, dis-assign the stereopermutator
    assignment_ = boost::none;
  }
}

void BondStereopermutator::Impl::propagateGraphChange(
  const AtomStereopermutator::PropagatedState& oldPermutatorState,
  const AtomStereopermutator& newPermutator,
  const PrivateGraph& graph,
  const StereopermutatorList& permutators
) {
  const RankingInformation& oldRanking = std::get<0>(oldPermutatorState);
  const AtomStereopermutator::ShapeMap& oldShapeMap = std::get<1>(oldPermutatorState);

  // We assume that the supplied permutators (or their state) were assigned
  assert(!oldShapeMap.empty());
  assert(newPermutator.assigned());

  /* We assume the old and new shapes are of the same size (i.e. this is
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

  // Recheck feasibility
  std::vector<unsigned> newFeasiblePermutations = notObviouslyInfeasibleStereopermutations(
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
   * - Since composite only offers dihedral information in terms of shape
   *   vertices, we have to translate these back into site indices
   * - Transform shape vertices through the sequence
   *
   *   old shape vertex
   *   -> site index (using old permutation state)
   *   -> atom index (using old ranking)
   *   -> site index (using new ranking)
   *   -> new shape vertex (using new permutation state)
   *
   *   and then compare dihedral values between the composites. This scheme
   *   has works even if the fused position has changed.
   * - Since newComposite may have reordered the OrientationStates, we have
   *   to be careful which part of the DihedralTuples we extract shape
   *   vertices from.
   * - Inversions of the dihedral defining sequence do not invert the sign of
   *   the dihedral
   */
  bool modifiedOrientationIsFirstInNewComposite = (
    newComposite.orientations().first.identifier
    == possiblyModifiedOrientation.identifier
  );

  using DihedralTuple = Stereopermutations::Composite::Permutation::DihedralTuple;
  // We know the permutator is assigned from a few lines earlier
  const std::vector<DihedralTuple>& oldDihedralList = composite_.allPermutations().at(
    feasiblePermutations_.at(*assignment_)
  ).dihedrals;

  auto getNewShapeVertex = [&](Shapes::Vertex oldVertex) -> Shapes::Vertex {
    const SiteIndex oldSiteIndex = oldShapeMap.indexOf(oldVertex);
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
   *
   * Watch out here: The exact ordering of shape vertices in the dihedral list
   * is not reproduced because we don't reorder angle vertices according to
   * their associated ranking character as in Composite. So we need set
   * membership tests instead of lexicographic equality.
   */
  const auto matchIter = Temple::find_if(
    newComposite,
    [&](const auto& newCompositePermutation) -> bool {
      return newCompositePermutation.close(newCompositeDihedrals);
    }
  );

  // There should always be a match
  if(matchIter == std::end(newComposite)) {
    throw std::logic_error("Bug: no match found in new composite.");
  }

  const unsigned permutationIndex = matchIter - std::begin(newComposite);
  const unsigned basePermutationIndex = newComposite.rankingEquivalentBase(permutationIndex);
  const auto feasibleIter = Temple::find(newFeasiblePermutations, basePermutationIndex);
  if(feasibleIter == std::end(newFeasiblePermutations)) {
    throw std::logic_error("Bug: no match for permutation in feasible");
  }

  // Overwrite class state
  assignment_ = feasibleIter - std::begin(newFeasiblePermutations);
  composite_ = std::move(newComposite);
  feasiblePermutations_ = std::move(newFeasiblePermutations);
}

void BondStereopermutator::Impl::propagateVertexRemoval(const AtomIndex removedIndex) {
  using Stereopermutations::Composite;
  const auto updateIndex = [&](const AtomIndex index) -> AtomIndex {
    if(index > removedIndex) {
      return index - 1;
    }

    if(index == removedIndex) {
      return PrivateGraph::removalPlaceholder;
    }

    return index;
  };

  // Atom indices are used in composite identifiers and the placement edge
  composite_.updateIdentifiers(updateIndex);
  edge_ = BondIndex(updateIndex(edge_.first), updateIndex(edge_.second));
}

/* Information */
BondStereopermutator::Alignment BondStereopermutator::Impl::alignment() const {
  return static_cast<BondStereopermutator::Alignment>(composite_.alignment());
}

boost::optional<unsigned> BondStereopermutator::Impl::assigned() const {
  return assignment_;
}

bool BondStereopermutator::Impl::hasSameCompositeOrientation(const BondStereopermutator::Impl& other) const {
  return composite_ == other.composite_;
}

boost::optional<unsigned> BondStereopermutator::Impl::indexOfPermutation() const {
  return Temple::Optionals::map(
    assignment_,
    Temple::Functor::at(feasiblePermutations_)
  );
}

unsigned BondStereopermutator::Impl::numAssignments() const {
  return feasiblePermutations_.size();
}

unsigned BondStereopermutator::Impl::numStereopermutations() const {
  return composite_.countNonEquivalentPermutations();
}

std::string BondStereopermutator::Impl::info() const {
  using namespace std::string_literals;

  std::string returnString;

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
