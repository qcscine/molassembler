/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/DistanceGeometry/SpatialModel.h"

#include "boost/numeric/interval.hpp"
#include "boost/range/iterator_range_core.hpp"
#include "boost/graph/graphviz.hpp"
#include "Utils/Typenames.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Constants.h"
#include "CyclicPolygons.h"
#include "RingDecomposerLib.h"

#include "shapes/PropertyCaching.h"
#include "stereopermutation/Composites.h"
#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/SequentialPairs.h"
#include "temple/Adaptors/CyclicFrame.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"
#include "temple/Optionals.h"
#include "temple/SetAlgorithms.h"
#include "temple/Stringify.h"
#include "temple/Stl17.h"

#include "molassembler/Cycles.h"
#include "molassembler/Detail/Cartesian.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/Graph/PrivateGraph.h"
#include "molassembler/Log.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/Modeling/ShapeInference.h"
#include "molassembler/Molecule/MolGraphWriter.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/Stereopermutators/FeasiblePermutations.h"

#include <fstream>
#include <Eigen/Dense>

namespace Scine {
namespace molassembler {

template<typename ... Inds>
inline auto orderedSequence(Inds ... inds) {
  std::array<AtomIndex, sizeof...(inds)> indices {{
    static_cast<AtomIndex>(inds)...
  }};

  if(indices.front() > indices.back()) {
    std::reverse(indices.begin(), indices.end());
  }

  return indices;
}

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

namespace distance_geometry {

// General availability of static constexpr members
constexpr double SpatialModel::bondRelativeVariance;
constexpr double SpatialModel::angleRelativeVariance;
constexpr double SpatialModel::dihedralAbsoluteVariance;

constexpr ValueBounds SpatialModel::angleClampBounds;
const ValueBounds SpatialModel::defaultDihedralBounds {std::nextafter(-M_PI, 0), M_PI};

SpatialModel::SpatialModel(
  const Molecule& molecule,
  const Configuration& configuration
) : _molecule(molecule) {
  /* This is overall a pretty complicated constructor since it encompasses the
   * entire conversion from a molecular graph into some model of the internal
   * coordinates of all connected atoms, determining which conformations are
   * accessible.
   *
   * The rough sequence of operations is:
   * - Make helper variables
   * - Set 1-2 bounds.
   * - Gather information on local geometries of all non-terminal atoms, using
   *   the existing stereopermutator data and supplanting it with random-assignment
   *   inferred stereopermutators on all other non-terminal atoms or double bonds.
   * - Set internal angles of all small flat cycles
   * - Set all remaining 1-3 bounds with additional tolerance if atoms involved
   *   in the angle are part of a small cycle
   * - Add BondStereopermutator 1-4 bound information
   *
   * The manner in which fixed positions are incorporated is that, whenever an
   * internal coordinate is modeled (i.e. bond lengths, angles or dihedrals)
   * between all-fixed atoms, this property is extracted from the supplied
   * positions, not generated from graph or stereopermutator information.
   */

  // Check invariants
  if(
    molecule.stereopermutators().hasZeroAssignmentStereopermutators()
    || molecule.stereopermutators().hasUnassignedStereopermutators()
  ) {
    throw std::logic_error("Failed precondition: molecule has zero-assignment or unassigned stereopermutators");
  }

  for(AtomIndex i : boost::make_iterator_range(_molecule.graph().atoms())) {
    if(_molecule.graph().degree(i) == 1) {
      continue;
    }

    if(!_molecule.stereopermutators().option(i)) {
      throw std::logic_error("Non-terminal atom is missing an atom stereopermutator");
    }
  }

  // Helper variables
  const Cycles& cycleData = molecule.graph().cycles();
  auto smallestCycleMap = makeSmallestCycleMap(cycleData);

  // Check constraints on static constants
  static_assert(
    0.0 < bondRelativeVariance && bondRelativeVariance < 0.1,
    "SpatialModel static constant bond relative variance must fulfill"
    "0 < x < 0.1"
  );
  static_assert(
    0.0 < angleRelativeVariance && angleRelativeVariance < 0.1,
    "SpatialModel static constant angle absolute variance must fulfill"
    "0 < x < 0.1"
  );

  // Add all information pertaining to fixed positions immediately
  temple::forEach(
    temple::adaptors::allPairs(configuration.fixedPositions),
    [&](const auto& indexPositionPairA, const auto& indexPositionPairB) {
      double spatialDistance = cartesian::distance(
        indexPositionPairA.second,
        indexPositionPairB.second
      ) * Utils::Constants::angstrom_per_bohr ;

      _constraints.emplace(
        orderedSequence(indexPositionPairA.first, indexPositionPairB.first),
        ValueBounds {spatialDistance, spatialDistance}
      );
    }
  );

  // Add any fixed positions to an unordered map for fast access:
  FixedPositionsMapType fixedAngstromPositions;
  for(const auto& fixedPositionPair : configuration.fixedPositions) {
    fixedAngstromPositions.emplace(
      fixedPositionPair.first,
      fixedPositionPair.second * Utils::Constants::angstrom_per_bohr
    );
  }

  // Set bond distances
  _modelBondDistances(fixedAngstromPositions, configuration.spatialModelLoosening);

  /* For all flat cycles for which the internal angles can be determined
   * exactly, add that information
   */
  _modelFlatCycles(fixedAngstromPositions, configuration.spatialModelLoosening);

  /* Create a list of AtomStereopermutators who must emit chiral constraints
   * even if they are achiral in order to satisfy refinement expectations of
   * BondStereopermutators.
   */
  std::unordered_set<AtomIndex> forceConstraintEmissionSet;
  for(const auto& bondStereopermutator : _molecule.stereopermutators().bondStereopermutators()) {
    for(AtomIndex placedAtomIndex : bondStereopermutator.edge()) {
      forceConstraintEmissionSet.insert(placedAtomIndex);
    }
  }

  // Get 1-3 information from AtomStereopermutators
  for(const auto& stereopermutator : _molecule.stereopermutators().atomStereopermutators()) {
    addAtomStereopermutatorInformation(
      stereopermutator,
      _molecule.graph().inner(),
      configuration.spatialModelLoosening,
      fixedAngstromPositions,
      forceConstraintEmissionSet.count(stereopermutator.centralIndex()) > 0
    );
  }

  // Get 1-4 information from BondStereopermutators
  for(const auto& bondStereopermutator : _molecule.stereopermutators().bondStereopermutators()) {
    addBondStereopermutatorInformation(
      bondStereopermutator,
      _molecule.stereopermutators().option(
        bondStereopermutator.edge().first
      ).value(),
      _molecule.stereopermutators().option(
        bondStereopermutator.edge().second
      ).value(),
      configuration.spatialModelLoosening,
      fixedAngstromPositions
    );
  }

  /* Model spiro centers */
  _modelSpirocenters(fixedAngstromPositions);

  /* Add default angles and dihedrals for all adjacent sequences (unpopulated
   * if e.g. the internal coordinates could not be modelled) to avoid using vdw
   * radius sums as minimum bounds for connected atom sequences.
   */
  _addDefaultAngles();
  _addDefaultDihedrals();
}

void SpatialModel::setBondBoundsIfEmpty(
  std::array<AtomIndex, 2> bondIndices,
  ValueBounds bounds
) {
  // Check the precondition
  assert(bondIndices.front() < bondIndices.back());

  _bondBounds.emplace(
    std::move(bondIndices),
    std::move(bounds)
  );
}

void SpatialModel::setAngleBoundsIfEmpty(
  std::array<AtomIndex, 3> angleIndices,
  ValueBounds bounds
) {
  // Check the precondition
  assert(angleIndices.front() < angleIndices.back());

  _angleBounds.emplace(
    std::move(angleIndices),
    clamp(std::move(bounds), angleClampBounds)
  );
}

void SpatialModel::setDihedralBoundsIfEmpty(
  std::array<AtomIndex, 4> dihedralIndices,
  ValueBounds bounds
) {
  // Check the precondition
  assert(dihedralIndices.front() < dihedralIndices.back());

  _dihedralBounds.emplace(
    std::move(dihedralIndices),
    std::move(bounds)
  );
}

void SpatialModel::addAtomStereopermutatorInformation(
  const AtomStereopermutator& permutator,
  const PrivateGraph& graph,
  const double looseningMultiplier,
  const std::unordered_map<AtomIndex, Utils::Position>& fixedAngstromPositions,
  const bool forceChiralConstraintEmission
) {
  const auto& feasiblePermutations = permutator.getFeasible();
  const auto& ranking = permutator.getRanking();
  const AtomIndex centerAtom = permutator.centralIndex();

  /* Angle information addition: Rough outline
   * - For each site
   *   - Set the distance to the center
   *     - For multi-atom sites, correct the distance with the cone angle
   *   - Add angles between site-constituting atoms using the cone angle
   *     (if there are multiple site-constituting atoms at all)
   * - For each pair of sites
   *   - The angle bounds between sites are initially set by the idealized
   *     shape position angle modified by the upper cone angles
   *   - For each pair of atoms between both sites
   *     - The site angle bounds variance is modified by cycle multipliers
   *       and loosening
   */

  /* Concerning fixed positions, if the center isn't fixed, then we do not
   * have to do anything differently at all.
   */
  std::vector<bool> siteFixed;
  bool centerFixed = fixedAngstromPositions.count(centerAtom) > 0;
  if(centerFixed) {
    // Map from site indices to whether the entire site is fixed or not
    siteFixed = temple::map(
      ranking.sites,
      [&fixedAngstromPositions](const auto& siteAtomList) -> bool {
        return temple::all_of(
          siteAtomList,
          [&fixedAngstromPositions](const AtomIndex siteConstitutingIndex) -> bool {
            return fixedAngstromPositions.count(siteConstitutingIndex) > 0;
          }
        );
      }
    );
  }

  /* Intra-site modelling / Between atoms within each site */
  for(unsigned long siteI = 0; siteI < feasiblePermutations.siteDistances.size(); ++siteI) {
    /* Set the distance to the center:
     * If no cone information is present, do not correct the distance to the
     * site using the cone angle
     *
     * NOTE: This is probably superfluous as non-eta bonds are modelled by the
     * SpatialModel beforehand. Try removing this block below and see if that
     * causes any issues.
     */
    if(!feasiblePermutations.coneAngles.at(siteI)) {
      for(const AtomIndex i : ranking.sites.at(siteI)) {
        setBondBoundsIfEmpty(
          orderedSequence(i, centerAtom),
          feasiblePermutations.siteDistances.at(siteI)
        );
      }

      // No further work has to be done for single-atom sites
      continue;
    }

    if(centerFixed && siteFixed.at(siteI)) {
      // All center to site constituting atom distances are fixed *exactly*
      for(const AtomIndex i : ranking.sites.at(siteI)) {
        double bondDistance = cartesian::distance(
          fixedAngstromPositions.at(i),
          fixedAngstromPositions.at(centerAtom)
        );

        setBondBoundsIfEmpty(
          orderedSequence(i, centerAtom),
          ValueBounds {bondDistance, bondDistance}
        );
      }

      // and so are the site constituting atom to site constituting atom angles
      temple::forEach(
        temple::adaptors::allPairs(ranking.sites.at(siteI)),
        [&](const AtomIndex i, const AtomIndex j) {
          const double angle = cartesian::angle(
            fixedAngstromPositions.at(i),
            fixedAngstromPositions.at(centerAtom),
            fixedAngstromPositions.at(j)
          );

          setAngleBoundsIfEmpty(
            orderedSequence(i, centerAtom, j),
            ValueBounds {angle, angle}
          );
        }
      );
    } else {
      /* Distance of every site site atom index to the central atom assumptions
       * - Every haptic index is on the cone base circle
       * - Cone height is defined by feasiblePermutations.siteDistance
       * - Cone angle is defined by feasiblePermutations.coneAngle
       */
      const ValueBounds& coneAngleBounds = feasiblePermutations.coneAngles.at(siteI).value();

      const ValueBounds hypotenuseBounds {
        feasiblePermutations.siteDistances.at(siteI).lower / std::cos(coneAngleBounds.upper),
        feasiblePermutations.siteDistances.at(siteI).upper / std::cos(coneAngleBounds.lower)
      };

      // Set bond distance for each site member to hypotenuse bounds
      for(const AtomIndex i : ranking.sites.at(siteI)) {
        setBondBoundsIfEmpty(
          orderedSequence(i, centerAtom),
          hypotenuseBounds
        );
      }

      /* Set angles between site-constituting atoms within a single site
       * - Minimally 0° (if there were a zero-length bond)
       *   You could compute shortest possible bond constexpr and insert a trig
       *   calc here, but the bond level distance is supplied elsewhere by
       *   SpatialModel anyway, no need to duplicate that information
       * - Maximally 2 * the upper cone angle (but not more than M_PI)
       */
      temple::forEach(
        temple::adaptors::allPairs(ranking.sites.at(siteI)),
        [&](const AtomIndex i, const AtomIndex j) {
          setAngleBoundsIfEmpty(
            orderedSequence(i, centerAtom, j),
            distance_geometry::ValueBounds {
              0,
              std::min(M_PI, 2 * coneAngleBounds.upper)
            }
          );
        }
      );
    }
  }

  /* Inter-site modelling / Between sites */
  /* If for either site no cone angles could be calculated (currently only
   * happens if a haptic site does not match the existing modeling patterns),
   * we have to skip this step entirely and hope that the remaining modeling
   * can pick up the slack.
   *
   * NOTE: Cone angles are calculated for non-haptic sites too -> (0, 0).
   */
  for(unsigned long i = 0; i < ranking.sites.size() - 1; ++i) {
    if(!feasiblePermutations.coneAngles.at(i)) {
      continue;
    }

    for(unsigned long j = i + 1; j < ranking.sites.size(); ++j) {
      if(!feasiblePermutations.coneAngles.at(j)) {
        continue;
      }

      if(centerFixed && siteFixed.at(i) && siteFixed.at(j)) {
        // All angles are known exactly!
        temple::forEach(
          temple::adaptors::allPairs(
            ranking.sites.at(i),
            ranking.sites.at(j)
          ),
          [&](const AtomIndex x, const AtomIndex y) -> void {
            const double angle = cartesian::angle(
              fixedAngstromPositions.at(x),
              fixedAngstromPositions.at(centerAtom),
              fixedAngstromPositions.at(y)
            );

            setAngleBoundsIfEmpty(
              orderedSequence(x, centerAtom, y),
              ValueBounds {angle, angle}
            );
          }
        );
      } else {
        /* The computed angle bounds are valid for each pair of atoms
         * constituting each site
         */
        temple::forEach(
          temple::adaptors::allPairs(
            ranking.sites.at(i),
            ranking.sites.at(j)
          ),
          [&](const AtomIndex x, const AtomIndex y) -> void {
            setAngleBoundsIfEmpty(
              orderedSequence(x, centerAtom, y),
              modelSiteAngleBounds(
                permutator,
                {i, j},
                looseningMultiplier,
                graph
              )
            );
          }
        );
      }
    }
  }

  // Add chiral constraints
  for(
    const AtomStereopermutator::MinimalChiralConstraint& minimalConstraint :
    permutator.minimalChiralConstraints(forceChiralConstraintEmission)
  ) {
    _chiralConstraints.emplace_back(
      makeChiralConstraint(minimalConstraint, permutator, looseningMultiplier)
    );
  }

  /* Add weak (low weight) planarity-enforcing chiral constraints if the shape
   * is planar and has more than two vertices
   */
  const shapes::Shape shape = permutator.getShape();
  const unsigned S = shapes::size(shape);
  if(!shapes::threeDimensional(shape) && S > 2) {
    constexpr double tolerance = 0.1;
    constexpr double weight = 0.01;

    auto siteIndices = temple::iota<unsigned>(S);

    const auto& sites = permutator.getRanking().sites;
    assert(sites.size() == S);
    for(unsigned offset = 0; offset < S - 2; ++offset) {
      ChiralConstraint::SiteSequence constraintSites;
      constraintSites[0] = {permutator.centralIndex()};
      constraintSites[1] = sites.at(offset);
      constraintSites[2] = sites.at(offset + 1);
      constraintSites[3] = sites.at(offset + 2);

      _chiralConstraints.emplace_back(
        std::move(constraintSites),
        -tolerance,
        tolerance
      );
      _chiralConstraints.back().weight = weight;
    }
  }
}

void SpatialModel::addBondStereopermutatorInformation(
  const BondStereopermutator& permutator,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const double looseningMultiplier,
  const std::unordered_map<AtomIndex, Utils::Position>& fixedAngstromPositions
) {
  // Check preconditions and get access to commonly needed things
  assert(permutator.indexOfPermutation());
  const unsigned permutation = permutator.indexOfPermutation().value();

  const stereopermutation::Composite& composite = permutator.composite();

  const AtomStereopermutator& firstStereopermutator = (
    stereopermutatorA.centralIndex() == composite.orientations().first.identifier
    ? stereopermutatorA
    : stereopermutatorB
  );

  const AtomStereopermutator& secondStereopermutator = (
    stereopermutatorB.centralIndex() == composite.orientations().second.identifier
    ? stereopermutatorB
    : stereopermutatorA
  );

  /* Only dihedrals */
  unsigned firstShapePosition;
  unsigned secondShapePosition;
  double dihedralAngle;

  for(const auto& dihedralTuple : composite.dihedrals(permutation)) {
    std::tie(firstShapePosition, secondShapePosition, dihedralAngle) = dihedralTuple;

    const unsigned siteIndexIAtFirst = SymmetryMapHelper::getSiteIndexAt(
      firstShapePosition,
      firstStereopermutator.getShapePositionMap()
    );
    const unsigned siteIndexLAtSecond = SymmetryMapHelper::getSiteIndexAt(
      secondShapePosition,
      secondStereopermutator.getShapePositionMap()
    );

    const auto& coneAngleIOption = firstStereopermutator.getFeasible().coneAngles.at(siteIndexIAtFirst);
    const auto& coneAngleLOption = secondStereopermutator.getFeasible().coneAngles.at(siteIndexLAtSecond);

    // Do not emit chiral constraints if cone angles are unknown
    if(!coneAngleIOption || !coneAngleLOption) {
      continue;
    }

    const ValueBounds& coneAngleI = *coneAngleIOption;
    const ValueBounds& coneAngleL = *coneAngleLOption;

    double dihedralVariance = coneAngleI.upper + coneAngleL.upper;
    if(permutator.alignment() == BondStereopermutator::Alignment::Eclipsed) {
      dihedralVariance += dihedralAbsoluteVariance * looseningMultiplier;
    } else if(permutator.alignment() == BondStereopermutator::Alignment::Staggered) {
      // Staggered dihedrals can be significantly looser
      dihedralVariance += 5 * dihedralAbsoluteVariance * looseningMultiplier;
    }

    /* If the width of the dihedral angle is now larger than 2π, then we may
     * overrepresent some dihedral values when choosing randomly in that
     * interval, and it is preferable just not to emit a dihedral constraint or
     * enter any dihedral distance information (the default values are covered
     * by addDefaultDihedrals).
     *
     * This should be very rare or not occur at all; it's just a safeguard.
     */
    if(dihedralVariance >= M_PI) {
      continue;
    }

    /* Modify the dihedral angle by the upper cone angles of the i and l
     * sites and the usual variances.
     *
     * NOTE: Don't worry about the periodicity issue here, the error function
     * terms in refinement have to take care of that.
     */
    ValueBounds dihedralBounds = makeBoundsFromCentralValue(
      dihedralAngle,
      dihedralVariance
    );

    // Set per-atom sequence dihedral distance bounds
    temple::forEach(
      temple::adaptors::allPairs(
        firstStereopermutator.getRanking().sites.at(siteIndexIAtFirst),
        secondStereopermutator.getRanking().sites.at(siteIndexLAtSecond)
      ),
      [&](const AtomIndex firstIndex, const AtomIndex secondIndex) -> void {
        setDihedralBoundsIfEmpty(
          orderedSequence(
            firstIndex,
            firstStereopermutator.centralIndex(),
            secondStereopermutator.centralIndex(),
            secondIndex
          ),
          dihedralBounds
        );
      }
    );

    /* Depending on alignment, we do not want to overdo the number of dihedral
     * constraints. If alignment is staggered, then we only do one-to-all
     * dihedrals instead of all-to-all.
     */
    if(
      composite.alignment() == stereopermutation::Composite::Alignment::Staggered
      && std::get<0>(composite.dihedrals(permutation).front()) != firstShapePosition
    ) {
      continue;
    }

    _dihedralConstraints.emplace_back(
      DihedralConstraint::SiteSequence {
        firstStereopermutator.getRanking().sites.at(siteIndexIAtFirst),
        {firstStereopermutator.centralIndex()},
        {secondStereopermutator.centralIndex()},
        secondStereopermutator.getRanking().sites.at(siteIndexLAtSecond)
      },
      dihedralBounds.lower,
      dihedralBounds.upper
    );
  }
}

template<std::size_t N>
bool bondInformationIsPresent(
  const DistanceBoundsMatrix& bounds,
  const std::array<AtomIndex, N>& indices
) {
  // Ensure all indices are unique
  std::set<AtomIndex> indicesSet {indices.begin(), indices.end()};
  if(indicesSet.size() < indices.size()) {
    return false;
  }

  // Check that the bond information in the sequence is present
  return temple::all_of(
    temple::adaptors::sequentialPairs(indices),
    [&bounds](const AtomIndex i, const AtomIndex j) -> bool {
      return (
        bounds.lowerBound(i, j) != DistanceBoundsMatrix::defaultLower
        && bounds.upperBound(i, j) != DistanceBoundsMatrix::defaultUpper
      );
    }
  );
}

SpatialModel::BoundsMatrix SpatialModel::makePairwiseBounds(
  unsigned N,
  const BoundsMapType<2>& fixedPositionBounds,
  const BoundsMapType<2>& bondBounds,
  const BoundsMapType<3>& angleBounds,
  const BoundsMapType<4>& dihedralBounds
) {
  BoundsMatrix bounds(N, N);
  bounds.setZero();

  /* There may be overlapping and possibly conflicting information present in
   * the gathered data. If data affects the same atom-pair, we must ensure that
   * we merely raise the lower bound and lower the upper bound, while never
   * inverting the bounds overall.
   */
  auto addInformation = [&bounds](
    AtomIndex i,
    AtomIndex j,
    const ValueBounds& newBounds
  ) {
    assert(i != j);
    // Ensure i < j
    if(j < i) {
      std::swap(i, j);
    }

    double& lowerBound = bounds(j, i);
    double& upperBound = bounds(i, j);

    assert(newBounds.lower <= newBounds.upper);
    assert(lowerBound <= upperBound);

    if(lowerBound != 0.0 && upperBound != 0.0) {
      if(
        newBounds.lower > lowerBound
        && newBounds.lower < upperBound
      ) {
        lowerBound = newBounds.lower;
      }

      // Try to lower the upper bound
      if(
        newBounds.upper < upperBound
        && newBounds.upper > lowerBound
      ) {
        upperBound = newBounds.upper;
      }
    } else {
      lowerBound = newBounds.lower;
      upperBound = newBounds.upper;
    }
  };

  auto getBondBounds = [&bounds](
    AtomIndex i,
    AtomIndex j
  ) -> ValueBounds {
    if(j < i) {
      std::swap(i, j);
    }

    return ValueBounds {
      bounds(j, i),
      bounds(i, j)
    };
  };

  auto copyMapIntoBounds = [&](const BoundsMapType<2>& boundsMap) {
    for(const auto& indexArrayBoundsPair : boundsMap) {
      const std::array<AtomIndex, 2>& indexArray = indexArrayBoundsPair.first;
      const ValueBounds& valueBounds = indexArrayBoundsPair.second;

      const AtomIndex i = indexArray.front();
      const AtomIndex j = indexArray.back();

      assert(i < j);
      bounds(j, i) = valueBounds.lower;
      bounds(i, j) = valueBounds.upper;
    }
  };

  // Copy the constraints as ground truth
  copyMapIntoBounds(fixedPositionBounds);

  // Add 1-2 information from the bonds
  if(fixedPositionBounds.empty()) {
    /* If there are no ground constraints, then we can just copy over the
     * entire bond bounds (these are definitely compatible with triangle
     * inequalities)
     */
    copyMapIntoBounds(bondBounds);
  } else {
    // Otherwise, we have to carefully add bond information
    for(const auto& bondPair : bondBounds) {
      addInformation(
        bondPair.first.front(),
        bondPair.first.back(),
        bondPair.second
      );
    }
  }

  // Add 1-3 information
  for(const auto& anglePair : angleBounds) {
    const auto& indices = anglePair.first;
    const auto& angleValueBounds = anglePair.second;

    ValueBounds firstBounds = getBondBounds(indices.front(), indices.at(1));
    ValueBounds secondBounds = getBondBounds(indices.at(1), indices.back());

    addInformation(
      indices.front(),
      indices.back(),
      ValueBounds {
        CommonTrig::lawOfCosines(
          firstBounds.lower,
          secondBounds.lower,
          angleValueBounds.lower
        ),
        CommonTrig::lawOfCosines(
          firstBounds.upper,
          secondBounds.upper,
          angleValueBounds.upper
        )
      }
    );
  }

  // Add 1-4 information
  for(const auto& dihedralPair : dihedralBounds) {
    const auto& indices = dihedralPair.first;
    const auto& dihedralValueBounds = dihedralPair.second;

    ValueBounds firstBounds = getBondBounds(indices.front(), indices.at(1));
    ValueBounds secondBounds = getBondBounds(indices.at(1), indices.at(2));
    ValueBounds thirdBounds = getBondBounds(indices.at(2), indices.back());

    auto firstAngleFindIter = angleBounds.find(
      orderedSequence(
        indices.at(0),
        indices.at(1),
        indices.at(2)
      )
    );

    auto secondAngleFindIter = angleBounds.find(
      orderedSequence(
        indices.at(1),
        indices.at(2),
        indices.at(3)
      )
    );

    if(
      firstAngleFindIter == angleBounds.end()
      || secondAngleFindIter == angleBounds.end()
    ) {
      continue;
    }

    const auto& abAngleBounds = firstAngleFindIter->second;
    const auto& bcAngleBounds = secondAngleFindIter->second;

    addInformation(
      indices.front(),
      indices.back(),
      CommonTrig::dihedralLengthBounds(
        firstBounds,
        secondBounds,
        thirdBounds,
        abAngleBounds,
        bcAngleBounds,
        dihedralValueBounds
      )
    );
  }

  return bounds;
}

double SpatialModel::siteCentralAngle(
  const AtomIndex centralIndex,
  const shapes::Shape& shape,
  const RankingInformation& ranking,
  const std::vector<unsigned>& shapeVertexMap,
  const std::pair<unsigned, unsigned>& sites,
  const PrivateGraph& inner
) {
  /* We need to respect the graph as ground truth and distort an ideal
   * shape angle towards lower angles for small cycles under specific
   * circumstances.
   */
  const double idealAngle = shapes::angleFunction(shape)(
    shapeVertexMap.at(sites.first),
    shapeVertexMap.at(sites.second)
  );

  /* The shape does not distort if:
   *
   * - Either site is haptic: Additional rotational degrees of freedom for the
   *   haptic sites ought to be sufficient to enable cycles.
   * - The ideal angle between the two sites is not the shape's smallest
   *   angle. Otherwise, heavily distorted trans-arrangements may be determined
   *   as feasible even though much more favorable cis-arrangements are
   *   possible.
   */
  if(
    ranking.sites.at(sites.first).size() > 1
    || ranking.sites.at(sites.second).size() > 1
    || idealAngle != shapes::minimumAngle(shape)
  ) {
    return idealAngle;
  }

  /* The shape also does not distort if the central index isn't part of a
   * small cycle (i.e. size < 6). So we look for cycles that contain the central
   * index and the specified two monoatomic sites.
   */
  // Find cycles that contain the central index and its two monoatomic sites
  std::vector<BondIndex> prospectiveCycleEdges {
    BondIndex {
      centralIndex,
      ranking.sites.at(sites.first).front()
    },
    BondIndex {
      centralIndex,
      ranking.sites.at(sites.second).front()
    }
  };

  auto cycleIteratorPair = inner.cycles().containing(prospectiveCycleEdges);

  // If the range is zero-length, there are no cycles with both edges!
  if(cycleIteratorPair.first == cycleIteratorPair.second) {
    return idealAngle;
  }

  unsigned smallestCycleSize = 100;
  std::vector<BondIndex> minimalCycle;
  for(auto cycleEdges : boost::make_iterator_range(cycleIteratorPair)) {
    if(cycleEdges.size() < smallestCycleSize) {
      minimalCycle = std::move(cycleEdges);
      smallestCycleSize = minimalCycle.size();
    }
  }

  // If the smallest cycle isn't small, then the angle doesn't distort
  if(smallestCycleSize >= 6) {
    return idealAngle;
  }

  /* Now it's time to model how the angle distorts for the given cycle.
   *
   * For cycles of size three, the angle is known exactly, since the triangle
   * is fully determined by the bond lengths
   */
  if(smallestCycleSize == 3) {
    BondIndex lastEdge {
      ranking.sites.at(sites.first).front(),
      ranking.sites.at(sites.second).front()
    };

    return CommonTrig::lawOfCosinesAngle(
      modelDistance(prospectiveCycleEdges.front(), inner),
      modelDistance(prospectiveCycleEdges.back(), inner),
      modelDistance(lastEdge, inner)
    );
  }

  /* For cycles of size four and five, the situation is more complex.
   *
   * Best would be modeling the full cycle with DG and use triangle / tetrangle
   * inequality smoothing to get at the distance between the site atoms and
   * then calculate the angle bounds and average.
   *
   * Next best, and a lot easier, could be the angle of a cyclic polygon model.
   */

  // Make sure the bonds adjacent to the permutator are at the front of minimalCycle
  auto iterForSwap = std::begin(minimalCycle);
  for(auto it = std::begin(minimalCycle); it != std::end(minimalCycle); ++it) {
    if(it->first == centralIndex || it->second == centralIndex) {
      std::iter_swap(iterForSwap, it);
      ++iterForSwap;
    }
  }

  // Model the cyclic polygon
  auto internalAngles = cyclic_polygons::internalAngles(
    temple::map(
      minimalCycle,
      [&](const BondIndex& bond) -> double {
        return modelDistance(bond, inner);
      }
    )
  );

  /* The first angle is that between the first two edge lengths (which we have
   * ensured are those edge lengths of the bonds adjacent to the permutator)
   */
  return internalAngles.front();
}

double SpatialModel::siteCentralAngle(
  const AtomStereopermutator& permutator,
  const std::pair<unsigned, unsigned>& sites,
  const PrivateGraph& inner
) {
  return siteCentralAngle(
    permutator.centralIndex(),
    permutator.getShape(),
    permutator.getRanking(),
    permutator.getShapePositionMap(),
    sites,
    inner
  );
}

ValueBounds SpatialModel::modelSiteAngleBounds(
  const AtomStereopermutator& permutator,
  const std::pair<unsigned, unsigned>& sites,
  const double looseningMultiplier,
  const PrivateGraph& inner
) {
  assert(permutator.assigned());

  const stereopermutators::Feasible& feasiblePermutations = permutator.getFeasible();

  /* The idealized shape angles are modified by the upper (!) cone angles
   * at each site, not split between lower and upper.
   */
  const double centralAngle = siteCentralAngle(permutator, sites, inner);
  const double absoluteVariance = [&]() -> double{
    // Turn the angle relative variance into an absolute variance
    double variance = SpatialModel::angleRelativeVariance * centralAngle;
    variance *= smallestCycleDistortionMultiplier(
      permutator.centralIndex(),
      inner.cycles()
    );
    variance *= looseningMultiplier;

    // Additional terms are the upper(!) cone angles!
    variance += feasiblePermutations.coneAngles.at(sites.first).value().upper;
    variance += feasiblePermutations.coneAngles.at(sites.second).value().upper;

    return variance;
  }();

  return clamp(
    makeBoundsFromCentralValue(centralAngle, absoluteVariance),
    angleClampBounds
  );
}

ChiralConstraint SpatialModel::makeChiralConstraint(
  const AtomStereopermutator::MinimalChiralConstraint& minimalConstraint,
  const AtomStereopermutator& permutator,
  const double looseningMultiplier
) {
  const auto& feasiblePermutations = permutator.getFeasible();
  const auto& ranking = permutator.getRanking();
  const AtomIndex centerAtom = permutator.centralIndex();

  /* We need to calculate target upper and lower volumes for the chirality
   * constraints. _cache.siteDistances contains bounds for the distance to
   * each site site plane, and since the center of each cone should
   * constitute the average site position, we can calculate 1-3 distances
   * between the centerpoints of sites using the idealized angles.
   *
   * The target volume of the chiral constraint created by the
   * tetrahedron is calculated using internal coordinates (the
   * Cayley-Menger determinant), always leading to V > 0, so depending on
   * the current assignment, the sign of the result is switched. The
   * formula used later in chiral constraint calculation for explicit
   * coordinates is adjusted by V' = 6 V to avoid an unnecessary factor, so
   * we do that here too:
   *
   *    288 V²  = |...|               | substitute V' = 6 V
   * -> 8 (V')² = |...|
   * ->      V' = sqrt(|...| / 8)
   *
   * where the Cayley-Menger determinant |...| is square symmetric:
   *
   *          |   0    1    1    1    1  |
   *          |        0  d12² d13² d14² |
   *  |...| = |             0  d23² d24² |
   *          |                  0  d34² |
   *          |  ...                  0  |
   *
   */

  constexpr unsigned cayleyMengerDim = 5;
  using DeterminantMatrix = Eigen::Matrix<double, cayleyMengerDim, cayleyMengerDim>;

  DeterminantMatrix lowerMatrix;
  DeterminantMatrix upperMatrix;

  lowerMatrix.row(0).setOnes();
  upperMatrix.row(0).setOnes();

  lowerMatrix.diagonal().setZero();
  upperMatrix.diagonal().setZero();

  /* Cycle through all combinations of site indices in the tetrahedron
   * definition sequence. boost::none means the central atom.
   */
  for(unsigned i = 0; i < 4; ++i) {
    boost::optional<distance_geometry::ValueBounds> iBounds;
    if(minimalConstraint.at(i)) {
      iBounds = feasiblePermutations.siteDistances.at(
        minimalConstraint.at(i).value()
      );
    }

    for(unsigned j = i + 1; j < 4; ++j) {
      boost::optional<distance_geometry::ValueBounds> jBounds;
      if(minimalConstraint.at(j)) {
        jBounds = feasiblePermutations.siteDistances.at(
          minimalConstraint.at(j).value()
        );
      }

      assert(iBounds || jBounds);

      distance_geometry::ValueBounds oneThreeDistanceBounds;
      if(iBounds && jBounds) {
        /* If neither index is the central atom, we can calculate an
         * expected one-three distance
         */
        double siteAngle = permutator.angle(
          minimalConstraint.at(i).value(),
          minimalConstraint.at(j).value()
        );

        oneThreeDistanceBounds = {
          CommonTrig::lawOfCosines(
            iBounds.value().lower,
            jBounds.value().lower,
            std::max(0.0, siteAngle * (1 - angleRelativeVariance * looseningMultiplier))
          ),
          CommonTrig::lawOfCosines(
            iBounds.value().upper,
            jBounds.value().upper,
            std::min(M_PI, siteAngle * (1 + angleRelativeVariance * looseningMultiplier))
          )
        };
      } else if(iBounds) {
        oneThreeDistanceBounds = iBounds.value();
      } else {
        oneThreeDistanceBounds = jBounds.value();
      }

      lowerMatrix(i + 1, j + 1) = std::pow(oneThreeDistanceBounds.lower, 2);
      upperMatrix(i + 1, j + 1) = std::pow(oneThreeDistanceBounds.upper, 2);
    }
  }

  const double boundFromLower = static_cast<DeterminantMatrix>(
    lowerMatrix.selfadjointView<Eigen::Upper>()
  ).determinant();

  const double boundFromUpper = static_cast<DeterminantMatrix>(
    upperMatrix.selfadjointView<Eigen::Upper>()
  ).determinant();

  assert(boundFromLower > 0 && boundFromUpper > 0);

  const double volumeFromLower = std::sqrt(boundFromLower / 8);
  const double volumeFromUpper = std::sqrt(boundFromUpper / 8);

  // Map the site indices to their constituent indices for use in the prototype
  auto tetrahedronSites = temple::map(
    minimalConstraint,
    [&](const boost::optional<unsigned>& siteIndexOptional) -> std::vector<AtomIndex> {
      if(siteIndexOptional) {
        return ranking.sites.at(siteIndexOptional.value());
      }

      return {centerAtom};
    }
  );

  /* Although it is tempting to assume that the Cayley-Menger determinant
   * using the lower bounds is smaller than the one using upper bounds,
   * this is not always true. We cannot a priori know which of both yields
   * the lower or upper bounds on the 3D volume, and hence must ensure only
   * that the ordering is preserved in the generation of the constraint,
   * which checks that the lower bound on the volume is smaller than the
   * upper one.
   *
   * You can check this assertion with a CAS. The relationship between both
   * determinants (where u_ij = l_ij + Δ) is wholly indeterminate, i.e. no
   * logical operator (<, >, <=, >=, ==) between both is true. It
   * completely depends on the individual values.
   *
   * Helpfully, chemical_symmetry only emits positive chiral target volume
   * index sequences (see test case name allTetrahedraPositive), so no
   * negative volumes have to be considered.
   */

  return {
    std::move(tetrahedronSites),
    std::min(volumeFromLower, volumeFromUpper),
    std::max(volumeFromLower, volumeFromUpper)
  };
}

SpatialModel::BoundsMatrix SpatialModel::makePairwiseBounds() const {
  return makePairwiseBounds(
    _molecule.graph().N(),
    _constraints,
    _bondBounds,
    _angleBounds,
    _dihedralBounds
  );
}

std::vector<distance_geometry::ChiralConstraint> SpatialModel::getChiralConstraints() const {
  return _chiralConstraints;
}

std::vector<distance_geometry::DihedralConstraint> SpatialModel::getDihedralConstraints() const {
  return _dihedralConstraints;
}

struct SpatialModel::ModelGraphWriter final : public MolGraphWriter {
  /* State */
  const SpatialModel& spatialModel;

/* Constructor */
  ModelGraphWriter(const PrivateGraph& inner, const SpatialModel& passSpatialModel)
    : MolGraphWriter(&inner, &passSpatialModel._molecule.stereopermutators()),
      spatialModel(passSpatialModel) {}

  std::vector<std::string> edgeTooltips(AtomIndex source, AtomIndex target) const final {
    const auto indexSequence = orderedSequence(source, target);
    if(spatialModel._bondBounds.count(indexSequence) == 1) {
      const auto& bondBounds = spatialModel._bondBounds.at(indexSequence);
      std::string tooltip = "[" + std::to_string(bondBounds.lower);
      tooltip += ", ";
      tooltip += std::to_string(bondBounds.upper);
      tooltip += "]";

      return {std::move(tooltip)};
    }

    return {};
  }

  std::vector<std::string> atomStereopermutatorTooltips(
    const AtomStereopermutator& permutator
  ) const final {
    std::vector<std::string> tooltips;
    tooltips.emplace_back(shapes::name(permutator.getShape()));
    tooltips.emplace_back(permutator.info());

    for(const auto& angleIterPair : spatialModel._angleBounds) {
      const auto& indexSequence = angleIterPair.first;
      const auto& angleBounds = angleIterPair.second;

      if(indexSequence.at(1) == permutator.centralIndex()) {
        tooltips.emplace_back(
          "["s + std::to_string(indexSequence.at(0)) + ","s
          + std::to_string(indexSequence.at(2)) +"] -> ["s
          + std::to_string(
            std::round(temple::Math::toDegrees(angleBounds.lower))
          ) + ", "s + std::to_string(
            std::round(temple::Math::toDegrees(angleBounds.upper))
          ) + "]"s
        );
      }
    }

    return tooltips;
  }

  std::vector<std::string> bondStereopermutatorTooltips(
    const BondStereopermutator& permutator
  ) const final {
    std::vector<std::string> tooltips;

    for(const auto& dihedralMapPair : spatialModel._dihedralBounds) {
      const auto& indexSequence = dihedralMapPair.first;
      const auto& dihedralBounds = dihedralMapPair.second;

      // Skip default dihedrals, list only explicitly set dihedrals
      if(dihedralBounds == defaultDihedralBounds) {
        continue;
      }

      if(
        (
          indexSequence.at(1) == permutator.edge().first
          && indexSequence.at(2) == permutator.edge().second
        ) || (
          indexSequence.at(1) == permutator.edge().second
          && indexSequence.at(2) == permutator.edge().first
        )
      ) {
        tooltips.emplace_back(
          "["s + std::to_string(indexSequence.at(0)) + ","s
          + std::to_string(indexSequence.at(3)) + "] -> ["s
          + std::to_string(
            std::round(temple::Math::toDegrees(dihedralBounds.lower))
          ) + ", "s + std::to_string(
            std::round(temple::Math::toDegrees(dihedralBounds.upper))
          ) + "]"s
        );
      }
    }

    return tooltips;
  }
};

std::string SpatialModel::dumpGraphviz() const {
  ModelGraphWriter graphWriter(
    _molecule.graph().inner(),
    *this
  );

  std::stringstream graphvizStream;

  boost::write_graphviz(
    graphvizStream,
    _molecule.graph().inner().bgl(),
    graphWriter,
    graphWriter,
    graphWriter
  );

  return graphvizStream.str();
}

void SpatialModel::writeGraphviz(const std::string& filename) const {
  ModelGraphWriter graphWriter(
    _molecule.graph().inner(),
    *this
  );

  std::ofstream outStream(filename);

  boost::write_graphviz(
    outStream,
    _molecule.graph().inner().bgl(),
    graphWriter,
    graphWriter,
    graphWriter
  );

  outStream.close();
}

/* Static functions */
double SpatialModel::modelDistance(
  const AtomIndex i,
  const AtomIndex j,
  const PrivateGraph& graph
) {
  return Bond::calculateBondDistance(
    graph.elementType(i),
    graph.elementType(j),
    graph.bondType(
      graph.edge(i, j)
    )
  );
}

double SpatialModel::modelDistance(
  const BondIndex& bond,
  const PrivateGraph& graph
) {
  return modelDistance(bond.first, bond.second, graph);
}

std::vector<BondIndex> SpatialModel::cycleConsistingOfExactly(
  const std::vector<AtomIndex>& atoms,
  const PrivateGraph& graph
) {
  std::vector<BondIndex> possibleCycleEdges;

  // Collect all graph edges between atoms of the cycle
  temple::forEach(
    temple::adaptors::allPairs(atoms),
    [&](const AtomIndex i, const AtomIndex j) {
      auto edgeOption = graph.edgeOption(i, j);
      if(edgeOption) {
        possibleCycleEdges.emplace_back(
          graph.source(*edgeOption),
          graph.target(*edgeOption)
        );
      }
    }
  );

  /* If the number of edges between the possible cycle atoms isn't the same
   * as the number of atoms, then it can't be a cycle anyway
   */
  if(possibleCycleEdges.size() != atoms.size()) {
    return {};
  }

  for(
    const auto cycleEdges :
    boost::make_iterator_range(
      graph.cycles().containing(possibleCycleEdges)
    )
  ) {
    if(cycleEdges.size() == atoms.size()) {
      return cycleEdges;
    }
  }

  return {};
}

boost::optional<ValueBounds> SpatialModel::coneAngle(
  const std::vector<AtomIndex>& baseConstituents,
  const ValueBounds& coneHeightBounds,
  const Graph& graph
) {
  /* Have to decide cone base radius in order to calculate this. There are some
   * simple cases to get out of the way first:
   */

  assert(!baseConstituents.empty());
  if(baseConstituents.size() == 1) {
    return ValueBounds {0.0, 0.0};
  }

  if(baseConstituents.size() == 2) {
    double radius = modelDistance(
      baseConstituents.front(),
      baseConstituents.back(),
      graph.inner()
    ) / 2;

    // Angle gets smaller if height bigger or cone base radius smaller
    double lowerAngle = std::atan2(
      (1 - bondRelativeVariance) * radius,
      coneHeightBounds.upper
    );

    double upperAngle = std::atan2(
      (1 + bondRelativeVariance) * radius,
      coneHeightBounds.lower
    );

    return ValueBounds {
      lowerAngle,
      upperAngle
    };
  }

  // Now it gets tricky. The base constituents may be part of a cycle or not
  auto cycleEdges = cycleConsistingOfExactly(baseConstituents, graph.inner());

  // The return value of cycleConsistingOfExactly is nullable
  if(!cycleEdges.empty()) {
    /* So if it IS a cycle, we need a ring index sequence to calculate a cyclic
     * polygon circumradius, which is how flat cycles are modelled here
     */
    auto ringIndexSequence = makeRingIndexSequence(
      std::move(cycleEdges)
    );

    auto distances = temple::map(
      temple::adaptors::cyclicFrame<2>(ringIndexSequence),
      [&](const AtomIndex i, const AtomIndex j) -> double {
        return modelDistance(i, j, graph.inner());
      }
    );

    auto lowerCircumradiusResult = cyclic_polygons::detail::convexCircumradius(
      temple::map(
        distances,
        [&](const double distance) -> double {
          return (1 - bondRelativeVariance) * distance;
        }
      )
    );

    auto upperCircumradiusResult = cyclic_polygons::detail::convexCircumradius(
      temple::map(
        distances,
        [&](const double distance) -> double {
          return (1 + bondRelativeVariance) * distance;
        }
      )
    );

    /* We assume that the circumcenter for any of these cyclic polygons should
     * be inside the polygon, not outside (meaning that the variation in edge
     * lengths is typically relatively small). If the circumcenter were outside
     * the polygon, it is clear that cycle atoms are not well approximated.
     */
    assert(upperCircumradiusResult.second);
    assert(lowerCircumradiusResult.second);

    return ValueBounds {
      std::atan2(lowerCircumradiusResult.first, coneHeightBounds.upper),
      std::atan2(upperCircumradiusResult.first, coneHeightBounds.lower)
    };
  }


  /* So the site atoms are NOT the sole constituents of a closed cycle.
   *
   * For some types of sites, we could still figure out a cone angle. If the
   * site group is actually a path in which any intermediate atom merely
   * connects the subsequent atoms (i.e. the path is not branched, there are no
   * cycles), and if all the involved intermediate geometries constituting the
   * longest path have only one distinct angle value, then we can create a
   * conformational model anyway.
   *
   * However, a path specific approach cannot treat branched haptic sites
   * (i.e. PN3 where both P and N bond to the metal), and we would need access
   * to the Molecule's StereopermutatorList in both cases. In that case, this
   * function, which should only be instrumental to devising which
   * stereopermutations are obviously impossible, is out of its depth. Perform
   * any additional modelling when the spatial model requires more information,
   * but not here.
   */
  return boost::none;
}

double SpatialModel::spiroCrossAngle(const double alpha, const double beta) {
  // The source of this equation is explained in documents/
  return std::acos(
    -std::cos(alpha / 2) * std::cos(beta / 2)
  );
}

ValueBounds SpatialModel::siteDistanceFromCenter(
  const std::vector<AtomIndex>& siteAtomList,
  const AtomIndex centralIndex,
  const Graph& graph
) {
  assert(!siteAtomList.empty());

  double distance;

  if(siteAtomList.size() == 1) {
    // Single-atom binding site
    AtomIndex atomIndex = siteAtomList.front();

    distance = modelDistance(
      atomIndex,
      centralIndex,
      graph.inner()
    );
  } else {
    // Haptic binding site
    distance = 0.9 * temple::average(
      temple::adaptors::transform(
        siteAtomList,
        [&](AtomIndex atomIndex) -> double {
          return modelDistance(
            atomIndex,
            centralIndex,
            graph.inner()
          );
        }
      )
    );
  }

  return {
    (1 - bondRelativeVariance) * distance,
    (1 + bondRelativeVariance) * distance
  };
}

double SpatialModel::smallestCycleDistortionMultiplier(
  const AtomIndex i,
  const Cycles& cycles
) {
  return temple::optionals::map(
    smallestCycleContaining(i, cycles),
    [](const unsigned cycleSize) -> double {
      if(cycleSize == 3) {
        return 6.25;
      }

      if(cycleSize == 4) {
        return 4.25;
      }

      if(cycleSize == 5) {
        return 3.25;
      }

      return 1.0;
    }
  ).value_or(1.0);
}

ValueBounds SpatialModel::makeBoundsFromCentralValue(
  const double centralValue,
  const double absoluteVariance
) {
  return {
    centralValue - absoluteVariance,
    centralValue + absoluteVariance
  };
}

ValueBounds SpatialModel::clamp(
  ValueBounds bounds,
  const ValueBounds& clampBounds
) {
  bounds.lower = temple::stl17::clamp(
    bounds.lower,
    clampBounds.lower,
    clampBounds.upper
  );

  bounds.upper = temple::stl17::clamp(
    bounds.upper,
    clampBounds.lower,
    clampBounds.upper
  );

  return bounds;
}

void SpatialModel::checkFixedPositionsPreconditions(
  const Molecule& molecule,
  const Configuration& configuration
) {
  // Early exit if there are no fixed positions
  if(configuration.fixedPositions.empty()) {
    return;
  }

  /* Check to ensure that every fixed atom has zero, one or all sites
   * fully fixed.
   */
  std::unordered_set<AtomIndex> fixedAtoms;
  for(const auto& fixedPositionPair : configuration.fixedPositions) {
    fixedAtoms.insert(fixedPositionPair.first);
  }

  for(const auto& indexPositionPair : configuration.fixedPositions) {
    const AtomIndex& atomIndex = indexPositionPair.first;

    auto stereopermutatorOption = molecule.stereopermutators().option(atomIndex);
    if(stereopermutatorOption) {
      // Check to ensure either 0, 1 or L sites are fixed
      unsigned numFixedSites = temple::accumulate(
        stereopermutatorOption->getRanking().sites,
        0U,
        [&fixedAtoms](const unsigned carry, const auto& indexSet) -> unsigned {
          const unsigned countFixed = temple::accumulate(
            indexSet,
            0U,
            [&fixedAtoms](const unsigned nestedCarry, const AtomIndex i) -> unsigned {
              return nestedCarry + fixedAtoms.count(i);
            }
          );

          if(countFixed != 0 && countFixed != indexSet.size()) {
            throw std::runtime_error(
              "DG preconditions for fixed atoms are not met: A non-terminal "
              "atom's binding site is only partially fixed."
            );
          }

          // Count this site if fully fixed
          if(countFixed == indexSet.size()) {
            return carry + 1;
          }

          return carry;
        }
      );

      if(1 < numFixedSites && numFixedSites < shapes::size(stereopermutatorOption->getShape())) {
        throw std::runtime_error(
          "DG preconditions for fixed atoms are not met: A non-terminal atom "
          "does not have 0, 1 or all binding sites fixed."
        );
      }
    }
  }
}

void SpatialModel::_addDefaultAngles() {
  const PrivateGraph& inner = _molecule.graph().inner();
  /* If no explicit angle can be provided for a triple of bonded atoms, we need
   * to at least specify the range of possible angles so that no implicit
   * minimimum distance (sum of vdw radii) is used instead. This is important
   * for haptic sites as for some site connectivity, no circumradius
   * can be modelled and hence angle calculation cannot be completed.
   */

  const AtomIndex N = _molecule.graph().N();
  for(AtomIndex center = 0; center < N; ++center) {
    temple::forEach(
      temple::adaptors::allPairs(
        boost::make_iterator_range(inner.adjacents(center))
      ),
      [&](const AtomIndex i, const AtomIndex j) -> void {
        assert(i != j);

        setAngleBoundsIfEmpty(
          orderedSequence(i, center, j),
          angleClampBounds
        );
      }
    );
  }
}

void SpatialModel::_addDefaultDihedrals() {
  const PrivateGraph& inner = _molecule.graph().inner();

  for(const auto& edgeDescriptor : boost::make_iterator_range(inner.edges())) {
    const AtomIndex sourceIndex = inner.source(edgeDescriptor);
    const AtomIndex targetIndex = inner.target(edgeDescriptor);

    temple::forEach(
      temple::adaptors::allPairs(
        boost::make_iterator_range(inner.adjacents(sourceIndex)),
        boost::make_iterator_range(inner.adjacents(targetIndex))
      ),
      [&](
        const AtomIndex sourceAdjacentIndex,
        const AtomIndex targetAdjacentIndex
      ) -> void {
        if(
          sourceAdjacentIndex != targetIndex
          && targetAdjacentIndex != sourceIndex
          && sourceAdjacentIndex != targetAdjacentIndex
        ) {
          setDihedralBoundsIfEmpty(
            orderedSequence(
              sourceAdjacentIndex,
              sourceIndex,
              targetIndex,
              targetAdjacentIndex
            ),
            defaultDihedralBounds
          );
        }
      }
    );
  }
}

void SpatialModel::_modelBondDistances(
  const FixedPositionsMapType& fixedAngstromPositions,
  const double looseningFactor
) {
  const PrivateGraph& inner = _molecule.graph().inner();

  for(const auto& edge: boost::make_iterator_range(inner.edges())) {
    BondType bondType = inner.bondType(edge);

    // Do not model eta bonds, stereopermutators are responsible for those
    if(bondType == BondType::Eta) {
      continue;
    }

    PrivateGraph::Vertex i = inner.source(edge);
    PrivateGraph::Vertex j = inner.target(edge);

    if(fixedAngstromPositions.count(i) > 0 && fixedAngstromPositions.count(j) > 0) {
      // If both atoms are fixed, their mutual bond distance is known exactly
      double bondDistance = cartesian::distance(
        fixedAngstromPositions.at(i),
        fixedAngstromPositions.at(j)
      );
      setBondBoundsIfEmpty(
        orderedSequence(i, j),
        ValueBounds {bondDistance, bondDistance}
      );
    } else {
      // Otherwise, create variable bounds on it
      double bondDistance = Bond::calculateBondDistance(
        inner.elementType(i),
        inner.elementType(j),
        bondType
      );

      double absoluteVariance = bondDistance * bondRelativeVariance * looseningFactor;

      setBondBoundsIfEmpty(
        orderedSequence(i, j),
        makeBoundsFromCentralValue(bondDistance, absoluteVariance)
      );
    }
  }
}

void SpatialModel::_modelFlatCycles(
  const FixedPositionsMapType& fixedAngstromPositions,
  const double looseningFactor
) {
  const PrivateGraph& inner = _molecule.graph().inner();
  const Cycles& cycleData = inner.cycles();

  for(auto cycleEdges : cycleData) {
    const unsigned cycleSize = cycleEdges.size();

    if(cycleSize >= 6) {
      continue;
    }

    /* There are a variety of cases here which all need to be treated
     * differently:
     * - Cycle of size 3: Always flat, set angles from cyclic polygons library
     * - Cycle of size 4: If flat (maybe already with a single double bond), use
     *   library to get precise internal angles and set them. Otherwise not all
     *   vertices are coplanar and all angles and dihedrals get a general
     *   tolerance increase
     * - Cycle of size 5: If aromatic, flat and we get precise internal angles
     *   from library. Otherwise, slightly increased angle and dihedral
     *   tolerances
     *
     * Assumptions
     * - Need to find out whether all vertices in cycles of size four are
     *   coplanar if one double bond is present. That seems to be the case for
     *   the few strained molecules collected in tests/strained_... Some more
     *   of these size four cycles are coplanar, especially when heteroatoms
     *   are present, but I don't immediately see the distinguishing criterion.
     */
    if( // flat cases
      cycleSize == 3
      || (
        cycleSize == 4
        && countPlanarityEnforcingBonds(
          cycleEdges,
          _molecule.graph()
        ) >= 1
      )
    ) {
      /* Gather sequence of atoms in cycle by progressively converting edge
       * descriptors into vertex indices
       */
      const auto indexSequence = makeRingIndexSequence(std::move(cycleEdges));

      // Skip any cycles that have fixed atoms
      if(
        temple::any_of(
          indexSequence,
          [&fixedAngstromPositions](const AtomIndex i) -> bool {
            return fixedAngstromPositions.count(i) > 0;
          }
        )
      ) {
        continue;
      }

      /* First, we fetch the angles that maximize the cycle area using the
       * cyclic polygon library.
       */
      const auto cycleInternalAngles = cyclic_polygons::internalAngles(
        // Map sequential index pairs to their purported bond length
        temple::map(
          temple::adaptors::cyclicFrame<2>(indexSequence),
          [&](const AtomIndex i, const AtomIndex j) -> double {
            return Bond::calculateBondDistance(
              inner.elementType(i),
              inner.elementType(j),
              inner.bondType(
                inner.edge(i, j)
              )
            );
          }
        )
      );

      /* The first angle returned is the angle between edges one and two, which
       * are indices [0, 1] and [1, 2] respectively. The last angle is between
       * edges [n - 1, 0], [0, 1].
       *
       * So, e.g. C4, C7, N9, O3::
       *
       *   index sequence 4, 7, 9, 3
       *   distances  C-C, C-N, N-O, O-C
       *   angles  C-C-N, C-N-O, N-O-C, O-C-C
       *                         ^- Last two overlap past index sequence
       *
       * Now we set all internal angles with low tolerance
       */
      // Make sure all assumptions we make about sequence lengths are true
      assert(indexSequence.size() == cycleInternalAngles.size());
      assert(indexSequence.size() == cycleSize);

      {
        /* This for loop construction is a little awkward, but it's really just
         * zipping cycleInternalAngles with the cyclic frame enumerating all
         * angle i-j-k sequences in the cycle.
         */
        unsigned angleIndex = 0;
        temple::forEach(
          temple::adaptors::cyclicFrame<3>(indexSequence),
          [&](const AtomIndex i, const AtomIndex j, const AtomIndex k) {
            setAngleBoundsIfEmpty(
              orderedSequence(i, j, k),
              makeBoundsFromCentralValue(
                cycleInternalAngles.at(angleIndex),
                angleRelativeVariance * looseningFactor * cycleInternalAngles.at(angleIndex)
              )
            );

            ++angleIndex;
          }
        );
      }
      /* Internal-external and external-external angles where the central
       * atom is part of a cycle are taken care of in the general angles
       * calculation below.
       */
    }
    /* Non-flat cases are also taken care of below using a general tolerance
     * increase on angles which directly effects the purported dihedral
     * distances.
     */
  }
}

void SpatialModel::_modelSpirocenters(
  const FixedPositionsMapType& fixedAngstromPositions
) {
  const PrivateGraph& inner = _molecule.graph().inner();
  const Cycles& cycleData = inner.cycles();

  /* For an atom to be a spiro center, it needs to be contained in exactly two
   * URFs and have four substituents in a tetrahedral symmetry.
   *
   * We also only have to worry about modeling them if the two URFs are
   * particularly small, i.e. are sizes 3-5
   */
  for(const auto& stereopermutator : _molecule.stereopermutators().atomStereopermutators()) {
    AtomIndex i = stereopermutator.centralIndex();

    // Skip any fixed central atoms
    if(fixedAngstromPositions.count(i) > 0) {
      continue;
    }

    // Skip any stereopermutators that do not match our conditions
    if(
      stereopermutator.getShape() != shapes::Shape::Tetrahedron
      || cycleData.numCycleFamilies(i) != 2
    ) {
      continue;
    }

    unsigned* URFIDs;
    auto nIDs = RDL_getURFsContainingNode(cycleData.dataPtr(), i, &URFIDs);
    assert(nIDs == 2);

    bool allURFsSingularRC = true;
    for(unsigned idIdx = 0; idIdx < nIDs; ++idIdx) {
      if(RDL_getNofRCForURF(cycleData.dataPtr(), URFIDs[idIdx]) > 1) {
        allURFsSingularRC = false;
        break;
      }
    }

    if(allURFsSingularRC) {
      // The two RCs still have to be disjoint save for i
      RDL_cycleIterator* cycleIteratorOne = RDL_getRCyclesForURFIterator(
        cycleData.dataPtr(),
        URFIDs[0]
      );
      RDL_cycle* cycleOne = RDL_cycleIteratorGetCycle(cycleIteratorOne);
      RDL_cycleIterator* cycleIteratorTwo = RDL_getRCyclesForURFIterator(
        cycleData.dataPtr(),
        URFIDs[1]
      );
      RDL_cycle* cycleTwo = RDL_cycleIteratorGetCycle(cycleIteratorTwo);

      // We model them only if both are small.
      if(cycleOne -> weight <= 5 && cycleTwo -> weight <= 5) {
        auto makeVerticesSet = [](const auto& cyclePtr) -> std::set<AtomIndex> {
          std::set<AtomIndex> vertices;

          for(unsigned cycleIndex = 0; cycleIndex < cyclePtr -> weight; ++cycleIndex) {
            vertices.insert(cyclePtr->edges[cycleIndex][0]);
            vertices.insert(cyclePtr->edges[cycleIndex][1]);
          }

          return vertices;
        };

        auto cycleOneVertices = makeVerticesSet(cycleOne);
        auto cycleTwoVertices = makeVerticesSet(cycleTwo);

        auto intersection = temple::set_intersection(
          cycleOneVertices,
          cycleTwoVertices
        );

        if(intersection.size() == 1 && *intersection.begin() == i) {
          std::vector<AtomIndex> firstAdjacents;
          std::vector<AtomIndex> secondAdjacents;
          for(
            const AtomIndex iAdjacent :
            boost::make_iterator_range(inner.adjacents(i))
          ) {
            if(cycleOneVertices.count(iAdjacent) > 0) {
              firstAdjacents.push_back(iAdjacent);
            }

            if(cycleTwoVertices.count(iAdjacent) > 0) {
              secondAdjacents.push_back(iAdjacent);
            }
          }

          assert(firstAdjacents.size() == 2 && secondAdjacents.size() == 2);

          auto firstSequence = orderedSequence(
            firstAdjacents.front(),
            i,
            firstAdjacents.back()
          );

          auto secondSequence = orderedSequence(
            secondAdjacents.front(),
            i,
            secondAdjacents.back()
          );

          /* We can only set these angles if the angle bounds for both
           * sequences already exist.
           */
          if(
            _angleBounds.count(firstSequence) == 1
            && _angleBounds.count(secondSequence) == 1
          ) {
            ValueBounds firstAngleBounds = _angleBounds.at(firstSequence);
            ValueBounds secondAngleBounds = _angleBounds.at(secondSequence);

            // Increases in cycle angles yield decrease in the cross angle
            double crossAngleLower = temple::stl17::clamp(
              spiroCrossAngle(
                firstAngleBounds.upper,
                secondAngleBounds.upper
              ),
              0.0,
              M_PI
            );

            double crossAngleUpper = temple::stl17::clamp(
              spiroCrossAngle(
                firstAngleBounds.lower,
                secondAngleBounds.lower
              ),
              0.0,
              M_PI
            );

            ValueBounds crossBounds {
              crossAngleLower,
              crossAngleUpper
            };

            temple::forEach(
              temple::adaptors::allPairs(
                firstAdjacents,
                secondAdjacents
              ),
              [&](const auto& firstAdjacent, const auto& secondAdjacent) {
                auto sequence = orderedSequence(firstAdjacent, i, secondAdjacent);

                auto findIter = _angleBounds.find(sequence);
                if(findIter == _angleBounds.end()) {
                  _angleBounds.emplace(
                    std::make_pair(sequence, crossBounds)
                  );
                } else {
                  // Tighten the cross-angle bounds
                  findIter->second = crossBounds;
                }
              }
            );
          }
        }
      }

      // Must manually free the cycles and iterators
      RDL_deleteCycle(cycleOne);
      RDL_deleteCycle(cycleTwo);

      RDL_deleteCycleIterator(cycleIteratorOne);
      RDL_deleteCycleIterator(cycleIteratorTwo);
    }

    // Must manually free the id array
    free(URFIDs);
  }
}

} // namespace distance_geometry
} // namespace molassembler
} // namespace Scine
