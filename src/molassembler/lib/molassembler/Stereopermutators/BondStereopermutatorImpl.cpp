/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/BondStereopermutatorImpl.h"

#include "shapes/Data.h"

#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/CyclicFrame.h"
#include "temple/constexpr/Math.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Functional.h"
#include "temple/OrderedPair.h"
#include "temple/Random.h"
#include "temple/STL17.h"

#include "molassembler/AngstromWrapper.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/Detail/Cartesian.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/ExplicitGraph.h"
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
namespace molassembler {
namespace detail {

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
 * of which the relevant types are std::pair and temple::OrderedPair here.
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

} // namespace detail

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

std::vector<char> BondStereopermutator::Impl::_charifyRankedSites(
  const RankingInformation::RankedSitesType& sitesRanking,
  const std::vector<unsigned>& symmetryPositionMap
) {
  std::vector<char> characters (symmetryPositionMap.size());

  char currentChar = 'A';
  for(const auto& equalPrioritySet : sitesRanking) {
    for(const auto& siteIndex : equalPrioritySet) {
      characters.at(
        symmetryPositionMap.at(siteIndex)
      ) = currentChar;
    }

    ++currentChar;
  }

  return characters;
}

const stereopermutation::Composite& BondStereopermutator::Impl::composite() const {
  return _composite;
}

double BondStereopermutator::Impl::dihedral(
  const AtomStereopermutator& stereopermutatorA,
  const unsigned siteIndexA,
  const AtomStereopermutator& stereopermutatorB,
  const unsigned siteIndexB
) const {
  if(!_assignment) {
    throw std::logic_error("This stereopermutator is unassigned! Dihedrals between sites are unspecified.");
  }

  bool swapped = false;

  // Use reference wrappers to be able to swap
  using PermutatorReference = std::reference_wrapper<const AtomStereopermutator>;
  using ReferencePair = std::pair<PermutatorReference, PermutatorReference>;

  ReferencePair references {stereopermutatorA, stereopermutatorB};
  std::pair<unsigned, unsigned> siteIndices {siteIndexA, siteIndexB};
  if(stereopermutatorA.centralIndex() == _composite.orientations().second.identifier) {
    std::swap(references.first, references.second);
    std::swap(siteIndices.first, siteIndices.second);
    swapped = true;
  }

  // Derived from possibly swapped references
  auto symmetryPositionMaps = temple::map_stl(
    references,
    [](const PermutatorReference& ref) {
      return ref.get().getShapePositionMap();
    }
  );

  std::pair<unsigned, unsigned> shapePositions;
  double dihedralAngle;

  for(const auto& dihedralTuple : _composite.dihedrals(*_assignment)) {
    std::tie(shapePositions.first, shapePositions.second, dihedralAngle) = dihedralTuple;

    if(
      SymmetryMapHelper::getSiteIndexAt(
        shapePositions.first,
        symmetryPositionMaps.first
      ) == siteIndices.first
      && SymmetryMapHelper::getSiteIndexAt(
        shapePositions.second,
        symmetryPositionMaps.second
      ) == siteIndices.second
    ) {
      return (swapped) ? -dihedralAngle : dihedralAngle;
    }
  }

  throw std::logic_error("Could not find a dihedral angle for the specified sites");
}

stereopermutation::Composite BondStereopermutator::Impl::_constructComposite(
  const StereopermutatorList& stereopermutators,
  const BondIndex edge,
  const Alignment alignment
) {
  assert(stereopermutators.option(edge.first) && stereopermutators.option(edge.second));

  const auto& stereopermutatorA = stereopermutators.option(edge.first).value();
  const auto& stereopermutatorB = stereopermutators.option(edge.second).value();

  return {
    _makeOrientationState(stereopermutatorA, stereopermutatorB),
    _makeOrientationState(stereopermutatorB, stereopermutatorA),
    static_cast<stereopermutation::Composite::Alignment>(alignment)
  };
}

stereopermutation::Composite::OrientationState
BondStereopermutator::Impl::_makeOrientationState(
  const AtomStereopermutator& focalStereopermutator,
  const AtomStereopermutator& attachedStereopermutator
) {
  return {
    focalStereopermutator.getShape(),
    SymmetryMapHelper::getSymmetryPositionOf(
      focalStereopermutator.getRanking().getSiteIndexOf(
        attachedStereopermutator.centralIndex()
      ),
      focalStereopermutator.getShapePositionMap()
    ),
    _charifyRankedSites(
      focalStereopermutator.getRanking().siteRanking,
      focalStereopermutator.getShapePositionMap()
    ),
    focalStereopermutator.centralIndex()
  };
}

bool BondStereopermutator::Impl::cycleObviouslyInfeasible(
  const InnerGraph& graph,
  const StereopermutatorList& stereopermutators,
  const AtomStereopermutator& firstStereopermutator,
  const AtomStereopermutator& secondStereopermutator,
  std::tuple<AtomIndex, AtomIndex, double> dihedral,
  const LinkInformation& link
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
  const AtomIndex j = firstStereopermutator.centralIndex();
  const AtomIndex k = secondStereopermutator.centralIndex();
  const AtomIndex l = std::get<1>(dihedral);
  const double phi = std::get<2>(dihedral);

  if(link.cycleSequence.size() == 3) {
    /* For triangles, modeling the cycle makes zero sense. As long as the
     * dihedral is approximately zero, the cycle must be feasible since the
     * graph is ground truth.
     */
    return std::fabs(phi) < temple::Math::toRadians(5.0);
  }

  if(link.cycleSequence.size() == 4) {
    /* It doesn't make much more sense to model the cycle for quadrangles, since
     * AD is actually bonded and extending that distance too much is essentially
     * breaking the bond. Make a little more allowance for the quadrangle
     * than the triangle, though.
     */
    return std::fabs(phi) < temple::Math::toRadians(20.0);
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

    const double lambda_m = temple::stl17::clamp(lambda, 0.0, 1.0);

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
  temple::inplace::remove_if(
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

  const auto elementTypes = temple::map(
    cycleIndices,
    [&](const AtomIndex x) -> Utils::ElementType {
      return graph.elementType(x);
    }
  );

  auto cycleEdgeLengths = temple::map(
    temple::adaptors::sequentialPairs(cycleIndices),
    modelDistance
  );
  cycleEdgeLengths.push_back((D-A).norm());

  // std::cout << "Cycle " << temple::stringify(link.cycleSequence)
  //   << "(a, b, c) = (" << a << ", " << b << ", " << c << "), "
  //   << "(alpha, beta, phi) = (" << temple::Math::toDegrees(alpha) << ", " << temple::Math::toDegrees(beta) << ", " << temple::Math::toDegrees(phi) << ") contradicts graph: " << fail << "\n";

  return Stereopermutators::cycleModelContradictsGraph(
    elementTypes,
    cycleEdgeLengths,
    bases
  );
}

// bool BondStereopermutator::Impl::cycleObviouslyInfeasible(
//   const InnerGraph& graph,
//   const StereopermutatorList& stereopermutators,
//   const AtomStereopermutator& firstStereopermutator,
//   const AtomStereopermutator& secondStereopermutator,
//   std::tuple<AtomIndex, AtomIndex, double> dihedral,
//   const LinkInformation& link
// ) {
//   /* Now to decide if a stereopermutation is possible or not:
//    * - Build a spatial model of each cycle including bond distances and angles
//    *   using SpatialModel's methods
//    * - Extract a BoundsList from it
//    * - Populate an ExplicitGraph
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
//   InnerGraph minimalInner(C);
//   // Copy element types
//   for(const AtomIndex i : link.cycleSequence) {
//     minimalInner.elementType(
//       indexReductionMap.at(i)
//     ) = graph.elementType(i);
//   }
//
//   // Copy bonds
//   temple::forEach(
//     temple::adaptors::cyclicFrame<2>(link.cycleSequence),
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
//   temple::forEach(
//     temple::adaptors::cyclicFrame<2>(link.cycleSequence),
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
//         detail::orderMappedSequence(indexReductionMap, i, j),
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
//   temple::forEach(
//     temple::adaptors::cyclicFrame<3>(link.cycleSequence),
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
//         detail::orderMappedSequence(indexReductionMap, i, j, k),
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
//   const AtomIndex j = firstStereopermutator.centralIndex();
//   const AtomIndex k = secondStereopermutator.centralIndex();
//   const AtomIndex l = std::get<1>(dihedral);
//
//   assert(temple::makeContainsPredicate(link.cycleSequence)(i));
//   assert(temple::makeContainsPredicate(link.cycleSequence)(l));
//
//   dihedralBounds.emplace(
//     detail::orderMappedSequence(indexReductionMap, i, j, k, l),
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
//   DistanceGeometry::ExplicitGraph boundsGraph {
//     minimalInner,
//     pairwiseBounds
//   };
//
//   // Try to smooth the distance bounds with triangle inequalities
//   auto distanceBoundsResult = boundsGraph.makeDistanceBounds();
//
//   /*if(distanceBoundsResult) {
//     std::cout << "Link " << temple::stringify(link.cycleSequence) << " is viable on " << j << ", " << k << " with dihedral " << dihedralAngle << "\n";
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
  const InnerGraph& graph,
  const StereopermutatorList& stereopermutators,
  const stereopermutation::Composite& composite
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
    return temple::iota<unsigned>(compositePermutations);
  }

  auto getDihedralInformation = [&](
    const std::vector<stereopermutation::Composite::DihedralTuple>& dihedrals,
    const LinkInformation& link
  ) -> boost::optional<std::tuple<AtomIndex, AtomIndex, double>> {
    unsigned firstSymmetryPosition;
    unsigned secondSymmetryPosition;
    double dihedralAngle;

    for(const auto& dihedralTuple : dihedrals) {
      std::tie(firstSymmetryPosition, secondSymmetryPosition, dihedralAngle) = dihedralTuple;

      const unsigned siteIndexIAtFirst = SymmetryMapHelper::getSiteIndexAt(
        firstSymmetryPosition,
        permutatorReferences.first.getShapePositionMap()
      );
      if(siteIndexIAtFirst != link.indexPair.first) {
        continue;
      }

      const unsigned siteIndexLAtSecond = SymmetryMapHelper::getSiteIndexAt(
        secondSymmetryPosition,
        permutatorReferences.second.getShapePositionMap()
      );
      if(siteIndexLAtSecond != link.indexPair.second) {
        continue;
      }

      const auto& iSite = permutatorReferences.first.getRanking().sites.at(siteIndexIAtFirst);
      const auto& lSite = permutatorReferences.second.getRanking().sites.at(siteIndexLAtSecond);
      if(iSite.size() > 1 || lSite.size() > 1) {
        return boost::none;
      }

      return std::make_tuple(
        iSite.front(),
        lSite.front(),
        dihedralAngle
      );
    }

    return boost::none;
  };

  std::vector<unsigned> viableStereopermutations;
  for(
    unsigned stereopermutationIndex = 0;
    stereopermutationIndex < compositePermutations;
    ++stereopermutationIndex
  ) {
    if(
      !temple::any_of(
        links,
        [&](const LinkInformation& link) -> bool {
          auto dihedralInformationOption = getDihedralInformation(
            composite.dihedrals(stereopermutationIndex),
            link
          );

          if(!dihedralInformationOption) {
            std::cout << "Could not match dihedral information\n";
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
) : _composite {
      _makeOrientationState(stereopermutatorA, stereopermutatorB),
      _makeOrientationState(stereopermutatorB, stereopermutatorA),
      static_cast<stereopermutation::Composite::Alignment>(alignment)
    },
    _edge(edge),
    _feasiblePermutations(
      temple::iota<unsigned>(_composite.permutations())
    ),
    _assignment(boost::none)
{}

BondStereopermutator::Impl::Impl(
  const InnerGraph& graph,
  const StereopermutatorList& stereopermutators,
  const BondIndex edge,
  Alignment alignment
) : _composite(
      _constructComposite(stereopermutators, edge, alignment)
    ),
    _edge(edge),
    _feasiblePermutations(
      notObviouslyInfeasibleStereopermutations(
        graph, stereopermutators, _composite
      )
    ),
    _assignment(boost::none)
{}

/* Public members */
/* Modification */
void BondStereopermutator::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment && assignment.value() >= numAssignments()) {
    /* The distinction here between numAssignments and _composite.permutations()
     * is important because if the composite is isotropic, we simulate that
     * there is only a singular assignment (see numAssignments), although
     * all permutations are generated and present for fitting anyway.
     *
     * If this were _composite.permutations(), we would accept assignments other
     * than zero if the underlying composite is isotropic, but not yield the
     * same assignment index when asked which assignment is currently set
     * in assigned().
     */
    throw std::out_of_range("Supplied assignment index is out of range");
  }

  _assignment = assignment;
}

void BondStereopermutator::Impl::assignRandom(random::Engine& engine) {
  const unsigned A = numAssignments();
  if(A == 0) {
    throw std::logic_error("Cannot randomly assign a stereopermutator without feasible stereopermutations");
  }

  if(A == 1) {
    assign(0);
  } else {
    assign(
      temple::random::getSingle<unsigned>(0, A - 1, engine)
    );
  }
}

void BondStereopermutator::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
  /* Composite's OrientationState identifiers change (these have no impact on
   * assignment ordering however, so that should be safe)
   */
  _composite.applyIdentifierPermutation(permutation);

  // Edge we sit on changes according to the permutation
  _edge = BondIndex {
    permutation.at(_edge.first),
    permutation.at(_edge.second)
  };

  // Assignment does not change
}

void BondStereopermutator::Impl::fit(
  const AngstromWrapper& angstromWrapper,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const FittingMode mode
) {
  assert(stereopermutatorA.centralIndex() != stereopermutatorB.centralIndex());

  // Early exit
  if(_composite.permutations() == 0) {
    _assignment = boost::none;
    return;
  }

  // Can the following selections be done with a single branch?
  const AtomStereopermutator& firstStereopermutator = (
    stereopermutatorA.centralIndex() == _composite.orientations().first.identifier
    ? stereopermutatorA
    : stereopermutatorB
  );

  const AtomStereopermutator& secondStereopermutator = (
    stereopermutatorB.centralIndex() == _composite.orientations().second.identifier
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

  unsigned firstSymmetryPosition;
  unsigned secondSymmetryPosition;
  double dihedralAngle;

  double bestPenalty = std::numeric_limits<double>::max();
  std::vector<unsigned> bestStereopermutation;

  for(unsigned feasiblePermutationIndex : _feasiblePermutations) {
    double penalty = 0.0;
    for(const auto& dihedralTuple : _composite.dihedrals(feasiblePermutationIndex)) {
      std::tie(firstSymmetryPosition, secondSymmetryPosition, dihedralAngle) = dihedralTuple;

      // Get site index of leftSymmetryPosition in left
      const unsigned firstSiteIndex = SymmetryMapHelper::getSiteIndexAt(
        firstSymmetryPosition,
        firstStereopermutator.getShapePositionMap()
      );
      const unsigned secondSiteIndex = SymmetryMapHelper::getSiteIndexAt(
        secondSymmetryPosition,
        secondStereopermutator.getShapePositionMap()
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

      const double measuredDihedral = cartesian::dihedral(
        firstSitePositions.col(firstSiteIndex),
        angstromWrapper.positions.row(firstStereopermutator.centralIndex()),
        angstromWrapper.positions.row(secondStereopermutator.centralIndex()),
        secondSitePositions.col(secondSiteIndex)
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
        penalty / _composite.dihedrals(feasiblePermutationIndex).size()
        > assignmentAcceptanceParameter * 2 * M_PI / _composite.order()
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
    _assignment = (
      temple::find(_feasiblePermutations, stereopermutation)
      - std::begin(_feasiblePermutations)
    );
  } else {
    // On ambiguous matching, dis-assign the stereopermutator
    _assignment = boost::none;
  }
}

void BondStereopermutator::Impl::propagateGraphChange(
  const AtomStereopermutatorPropagatedState& oldPermutatorState,
  const AtomStereopermutator& newPermutator,
  const InnerGraph& graph,
  const StereopermutatorList& permutators
) {
  const RankingInformation& oldRanking = std::get<0>(oldPermutatorState);
  const AbstractStereopermutations& oldAbstract = std::get<1>(oldPermutatorState);
  const FeasibleStereopermutations& oldFeasible = std::get<2>(oldPermutatorState);
  const boost::optional<unsigned>& oldAssignment = std::get<3>(oldPermutatorState);

  // We assume that the supplied permutators (or their state) were assigned
  assert(oldAssignment);
  assert(newPermutator.assigned());

  /* We assume the old and new symmetry are of the same size (i.e. this is
   * a ranking change propagation, not a substituent addition / removal
   * propagation)
   */
  assert(oldRanking.sites.size() == newPermutator.getRanking().sites.size());

  using OrientationState = stereopermutation::Composite::OrientationState;

  /* Construct a new Composite with the new information */
  bool changedIsFirstInOldOrientations = (
    _composite.orientations().first.identifier == newPermutator.centralIndex()
  );

  const OrientationState& oldOrientation = detail::select(
    _composite.orientations(),
    changedIsFirstInOldOrientations
  );

  // Reuse the OrientationState of the "other" atom stereopermutator
  const OrientationState& unchangedOrientation = detail::select(
    _composite.orientations(),
    !changedIsFirstInOldOrientations
  );

  // Generate a new OrientationState for the modified stereopermutator
  stereopermutation::Composite::OrientationState possiblyModifiedOrientation {
    newPermutator.getShape(),
    SymmetryMapHelper::getSymmetryPositionOf(
      newPermutator.getRanking().getSiteIndexOf(
        unchangedOrientation.identifier
      ),
      newPermutator.getShapePositionMap()
    ),
    _charifyRankedSites(
      newPermutator.getRanking().siteRanking,
      newPermutator.getShapePositionMap()
    ),
    newPermutator.centralIndex()
  };

  // In case nothing has changed, we are done and can stop
  if(oldOrientation == possiblyModifiedOrientation) {
    return;
  }

  /* Generate a new set of permutations (note this may reorder its
   * arguments into the orientations() member)
   */
  stereopermutation::Composite newComposite {
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
    _composite = std::move(newComposite);
    _feasiblePermutations = std::move(newFeasiblePermutations);
    if(!_feasiblePermutations.empty()) {
      _assignment = 0;
    }
    return;
  }

  /* If this BondStereopermutator is unassigned, and the permutator is not
   * isotropic after the ranking change, then this permutator stays unassigned.
   */
  if(_assignment == boost::none) {
    _composite = std::move(newComposite);
    _feasiblePermutations = std::move(newFeasiblePermutations);
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

  using DihedralTuple = stereopermutation::Composite::DihedralTuple;
  // This permutator is assigned since that is ensured a few lines earlier
  const std::vector<DihedralTuple>& oldDihedralList = _composite.dihedrals(
    _assignment.value()
  );

  auto oldSymmetryPositionToSiteMap = shapeVertexToSiteIndexMap(
    oldAbstract.permutations.stereopermutations.at(
      oldFeasible.indices.at(
        oldAssignment.value()
      )
    ),
    oldAbstract.canonicalSites
  );

  auto getNewSymmetryPosition = [&](unsigned oldSymmetryPosition) -> unsigned {
    const unsigned oldSiteIndex = oldSymmetryPositionToSiteMap.at(oldSymmetryPosition);

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

    const unsigned newSiteIndex = findSiteIter - std::begin(newRankingSites);

    return SymmetryMapHelper::getSymmetryPositionOf(
      newSiteIndex,
      newPermutator.getShapePositionMap()
    );
  };

  // Map the set of old dihedral tuples into the new space
  std::vector<DihedralTuple> newCompositeDihedrals;
  newCompositeDihedrals.reserve(oldDihedralList.size());
  for(const DihedralTuple& oldDihedral : oldDihedralList) {
    const unsigned changedSymmetryPosition = detail::select(
      oldDihedral,
      changedIsFirstInOldOrientations
    );

    const unsigned unchangedSymmetryPosition = detail::select(
      oldDihedral,
      !changedIsFirstInOldOrientations
    );

    const unsigned newSymmetryPosition = getNewSymmetryPosition(changedSymmetryPosition);

    if(modifiedOrientationIsFirstInNewComposite) {
      newCompositeDihedrals.emplace_back(
        newSymmetryPosition,
        unchangedSymmetryPosition,
        std::get<2>(oldDihedral)
      );
    } else {
      newCompositeDihedrals.emplace_back(
        unchangedSymmetryPosition,
        newSymmetryPosition,
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
            && detail::piPeriodicFPCompare(std::get<2>(a), std::get<2>(b))
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
  _assignment = matchIter - std::begin(newComposite);
  _composite = std::move(newComposite);
  _feasiblePermutations = std::move(newFeasiblePermutations);
}

/* Information */
BondStereopermutator::Alignment BondStereopermutator::Impl::alignment() const {
  return static_cast<BondStereopermutator::Alignment>(
    _composite.alignment()
  );
}

boost::optional<unsigned> BondStereopermutator::Impl::assigned() const {
  /* If the underlying composite is isotropic, it does not matter which of those
   * permutations by shape position is the factual spatial arrangement (since
   * they are all rotationally equivalent). We have to spoof that there is only
   * one arrangement in this case (although we need all of them for spatial
   * fitting).
   */
  if(_assignment && (_composite.isIsotropic() && !_feasiblePermutations.empty())) {
    return 0U;
  }

  return _assignment;
}

bool BondStereopermutator::Impl::hasSameCompositeOrientation(const BondStereopermutator::Impl& other) const {
  return _composite == other._composite;
}

boost::optional<unsigned> BondStereopermutator::Impl::indexOfPermutation() const {
  if(_assignment && _composite.isIsotropic()) {
    return 0U;
  }

  if(!_assignment) {
    return boost::none;
  }

  return _feasiblePermutations.at(*_assignment);
}

unsigned BondStereopermutator::Impl::numAssignments() const {
  if(_composite.isIsotropic() && !_feasiblePermutations.empty()) {
    return 1;
  }

  return _feasiblePermutations.size();
}

unsigned BondStereopermutator::Impl::numStereopermutations() const {
  if(_composite.isIsotropic()) {
    return 1;
  }

  return _composite.permutations();
}

std::string BondStereopermutator::Impl::info() const {
  using namespace std::string_literals;

  std::string returnString =  "B on ";

  returnString += std::to_string(_composite.orientations().first.identifier);
  returnString += "-";
  returnString += std::to_string(_composite.orientations().second.identifier);

  const unsigned A = numAssignments();

  if(A == 1) {
    returnString += ". Is non-stereogenic.";
  } else {
    returnString += ". Is ";
    if(_assignment) {
      returnString += std::to_string(_assignment.value());
    } else {
      returnString += "u";
    }

    returnString += "/"s + std::to_string(A);

    const unsigned P = numStereopermutations();
    if(P != A) {
      returnString += " ("s + std::to_string(P) + ")"s;
    }
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

BondIndex BondStereopermutator::Impl::edge() const {
  // Return a standard form of smaller first
  return _edge;
}

} // namespace molassembler

} // namespace Scine
