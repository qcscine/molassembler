// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/DistanceGeometry/SpatialModel.h"

#include "boost/numeric/interval.hpp"
#include "boost/range/iterator_range_core.hpp"
#include "boost/graph/graphviz.hpp"
#include "Delib/ElementInfo.h"
#include "Delib/Constants.h"
#include "CyclicPolygons.h"

#include "chemical_symmetries/Properties.h"
#include "stereopermutation/Composites.h"
#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/SequentialPairs.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"
#include "temple/Random.h"
#include "temple/SetAlgorithms.h"
#include "temple/Stringify.h"

#include "molassembler/Cycles.h"
#include "molassembler/Detail/DelibHelpers.h"
#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Log.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/Molecule/MolGraphWriter.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/Stereopermutators/PermutationState.h"

#include <fstream>
#include <Eigen/Dense>

namespace molassembler {

struct SymmetryMapHelper {
  static unsigned getSymmetryPositionOf(unsigned ligandIndex, const std::vector<unsigned>& map) {
    return map.at(ligandIndex);
  }

  static unsigned getLigandIndexAt(unsigned symmetryPosition, const std::vector<unsigned>& map) {
    auto findIter = std::find(
      std::begin(map),
      std::end(map),
      symmetryPosition
    );

    assert(findIter != std::end(map));

    return findIter - std::begin(map);
  }
};

namespace DistanceGeometry {

// General availability of static constexpr members
constexpr double SpatialModel::bondRelativeVariance;
constexpr double SpatialModel::angleAbsoluteVariance;
constexpr double SpatialModel::dihedralAbsoluteVariance;

constexpr ValueBounds SpatialModel::angleClampBounds;
ValueBounds SpatialModel::defaultDihedralBounds = {std::nextafter(-M_PI, 0), M_PI};

SpatialModel::SpatialModel(
  const Molecule& molecule,
  const Configuration& configuration,
  const double looseningMultiplier
) : _molecule(molecule),
    _looseningMultiplier(looseningMultiplier),
    _stereopermutators(molecule.stereopermutators())
{
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
   *
   *   - Copy original molecule's list of stereopermutators. Set unassigned ones to
   *     a random assignment (consistent with occurrence statistics)
   *   - Instantiate BondStereopermutators on all double bonds that aren't immediately
   *     involved in other stereopermutators. This is to ensure that all atoms
   *     involved in non-stereogenic double bonds are also flat in the
   *     resulting 3D structure. They additionally allow extraction of
   *     angle and dihedral angle information just like AtomStereopermutators.
   *   - Instantiate AtomStereopermutators on all remaining atoms and default-assign
   *     them. This is so that we can get angle data between substituents easily.
   *
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
  assert(
    !molecule.stereopermutators().hasZeroAssignmentStereopermutators()
    && !molecule.stereopermutators().hasUnassignedStereopermutators()
    && "The passed molecule may not have zero-assignment or unassigned stereopermutators!"
  );

  // Helper variables
  // Ignore eta bonds in construction of cycle data
  Cycles cycleData {molecule.graph(), true};
  auto smallestCycleMap = makeSmallestCycleMap(cycleData);

  // Check constraints on static constants
  static_assert(
    0 < bondRelativeVariance && bondRelativeVariance < 1,
    "SpatialModel static constant bond relative variance must fulfill"
    "0 < x << 1"
  );
  static_assert(
    0 < angleAbsoluteVariance && angleAbsoluteVariance < Symmetry::smallestAngle,
    "SpatialModel static constant angle absolute variance must fulfill"
    "0 < x << (smallest angle any local symmetry returns)"
  );

  const InnerGraph& inner = molecule.graph().inner();

  // Add all information pertaining to fixed positions immediately
  temple::forEach(
    temple::adaptors::allPairs(configuration.fixedPositions),
    [&](const auto& indexPositionPairA, const auto& indexPositionPairB) -> void {
      double spatialDistance = DelibHelpers::distance(
        indexPositionPairA.second.asEigenVector(),
        indexPositionPairB.second.asEigenVector()
      ) * Delib::angstrom_per_bohr ;

      _constraints.emplace(
        orderedIndexSequence<2>({{indexPositionPairA.first, indexPositionPairB.first}}),
        ValueBounds {spatialDistance, spatialDistance}
      );
    }
  );
  // Add any fixed positions to an unordered map for fast access:
  std::unordered_map<AtomIndex, Delib::Position> fixedAngstromPositions;

  for(const auto& fixedPositionPair : configuration.fixedPositions) {
    fixedAngstromPositions.emplace(
      fixedPositionPair.first,
      fixedPositionPair.second * Delib::angstrom_per_bohr
    );
  }

  // Set bond distances
  for(const auto& edge: boost::make_iterator_range(inner.edges())) {
    BondType bondType = inner.bondType(edge);

    // Do not model eta bonds, stereopermutators are responsible for those
    if(bondType == BondType::Eta) {
      continue;
    }

    InnerGraph::Vertex i = inner.source(edge);
    InnerGraph::Vertex j = inner.target(edge);

    if(fixedAngstromPositions.count(i) > 0 && fixedAngstromPositions.count(j) > 0) {
      // If both atoms are fixed, their mutual bond distance is known exactly
      double bondDistance = DelibHelpers::distance(
        fixedAngstromPositions.at(i).asEigenVector(),
        fixedAngstromPositions.at(j).asEigenVector()
      );
      setBondBoundsIfEmpty(
        {{i, j}},
        ValueBounds {bondDistance, bondDistance}
      );
    } else {
      // Otherwise, create variable bounds on it
      double bondDistance = Bond::calculateBondDistance(
        inner.elementType(i),
        inner.elementType(j),
        bondType
      );
      setBondBoundsIfEmpty(
        {{i, j}},
        bondDistance
      );
    }
  }

  /* The StereopermutatorList is already copy initialized with the Molecule's
   * stereopermutators, but we need to instantiate AtomStereopermutators everywhere,
   * regardless of whether they are stereogenic or not, to ensure that
   * modelling gets the information it needs.
   *
   * So for every missing non-terminal atom, create a AtomStereopermutator in the
   * determined geometry
   */
  const unsigned N = molecule.graph().N();
  for(unsigned i = 0; i < N; ++i) {
    // Already an instantiated AtomStereopermutator?
    if(_stereopermutators.option(i)) {
      continue;
    }

    auto localRanking = molecule.rankPriority(i);

    // Terminal atom index? skip those
    if(localRanking.ligands.size() <= 1) {
      continue;
    }

    Symmetry::Name localSymmetry = molecule.determineLocalGeometry(i, localRanking);

    auto newStereopermutator = AtomStereopermutator {
      molecule.graph(),
      localSymmetry,
      i,
      std::move(localRanking)
    };

    /* New stereopermutators encountered at this point can have multiple
     * assignments, since some types of stereopermutators are flatly ignored by the
     * candidate functions from Molecule, such as trigonal pyramidal nitrogens.
     * These are found here, though, and MUST be chosen randomly according to
     * the relative weights to get a single conformation in the final model
     */
    newStereopermutator.assignRandom();

    // Add it to the list of stereopermutators
    _stereopermutators.add(
      i,
      std::move(newStereopermutator)
    );
  }

  /* For all flat cycles for which the internal angles can be determined
   * exactly, add that information
   */
  for(
    const auto cyclePtr :
    boost::make_iterator_range(
      cycleData.iteratorPair(Cycles::predicates::SizeLessThan {6})
    )
  ) {
    const auto edgeDescriptors = Cycles::edges(cyclePtr);
    const unsigned cycleSize = edgeDescriptors.size();

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
          edgeDescriptors,
          molecule.graph()
        ) >= 1
      )
      /*! @todo missing cases:
       * - Need aromaticity checking routine for cycle size 5
       */
    ) {
      /* Gather sequence of atoms in cycle by progressively converting edge
       * descriptors into vertex indices
       */
      const auto indexSequence = makeRingIndexSequence(
        Cycles::edgeVertices(cyclePtr)
      );

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
      const auto cycleInternalAngles = CyclicPolygons::internalAngles(
        // Map sequential index pairs to their purported bond length
        temple::map(
          temple::adaptors::sequentialPairs(indexSequence),
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
       *   index sequence 4, 7, 9, 3, 4
       *   distances  C-C, C-N, N-O, O-C
       *   angles  C-C-N, C-N-O, N-O-C, O-C-C
       *                                ^- Last overlaps past index sequence
       *
       * Now we set all internal angles with low tolerance
       */
      // Make sure all assumptions we make about sequence lengths are true
      assert(indexSequence.size() == cycleInternalAngles.size() + 1);
      assert(indexSequence.size() == cycleSize + 1);
      // All non-overlapping triples
      for(
        unsigned angleCentralIndex = 1;
        angleCentralIndex < indexSequence.size() - 1;
        ++angleCentralIndex
      ) {
        setAngleBoundsIfEmpty(
          {{
            indexSequence.at(angleCentralIndex - 1),
            indexSequence.at(angleCentralIndex),
            indexSequence.at(angleCentralIndex + 1)
          }},
          makeBoundsFromCentralValue(
            cycleInternalAngles.at(angleCentralIndex - 1),
            angleAbsoluteVariance * looseningMultiplier
          )
        );
      }

      // There's a missing triple with the previous approach, which is always
      setAngleBoundsIfEmpty(
        {{
          indexSequence.at(indexSequence.size() - 2),
          indexSequence.at(0),
          indexSequence.at(1)
        }},
        makeBoundsFromCentralValue(
          cycleInternalAngles.back(),
          angleAbsoluteVariance * looseningMultiplier
        )
      );

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

  /* Returns a multiplier intended for the absolute angle variance for an atom
   * index. If that index is in a cycle of size < 6, the multiplier is > 1.
   */
  auto cycleMultiplierForIndex = [&](const AtomIndex i) -> double {
    auto findIter = smallestCycleMap.find(i);

    if(findIter != smallestCycleMap.end()) {
      if(findIter->second == 3) {
        return 6.25 ;
      }

      if(findIter->second == 4) {
        return 4.25;
      }

      if(findIter->second == 5) {
        return 3.25;
      }
    }

    return 1.0;
  };

  // Get 1-3 information from AtomStereopermutators
  for(const auto& stereopermutator : _stereopermutators.atomStereopermutators()) {
    addAtomStereopermutatorInformation(
      stereopermutator,
      cycleMultiplierForIndex,
      looseningMultiplier,
      fixedAngstromPositions
    );
  }

  // Get 1-4 information from BondStereopermutators
  for(const auto& bondStereopermutator : _stereopermutators.bondStereopermutators()) {
    addBondStereopermutatorInformation(
      bondStereopermutator,
      _stereopermutators.option(
        bondStereopermutator.edge().first
      ).value(),
      _stereopermutators.option(
        bondStereopermutator.edge().second
      ).value(),
      looseningMultiplier,
      fixedAngstromPositions
    );
  }

  /* Model spiro centers */
  /* For an atom to be a spiro center, it needs to be contained in exactly two
   * URFs and have four substituents in a tetrahedral symmetry.
   *
   * We also only have to worry about modeling them if the two URFs are
   * particularly small, i.e. are sizes 3-5
   */
  for(const auto& stereopermutator : _stereopermutators.atomStereopermutators()) {
    AtomIndex i = stereopermutator.centralIndex();

    // Skip any fixed central atoms
    if(fixedAngstromPositions.count(i) > 0) {
      continue;
    }

    // Skip any stereopermutators that do not match our conditions
    if(
      stereopermutator.getSymmetry() != Symmetry::Name::Tetrahedral
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
          std::vector<AtomIndex> firstAdjacents, secondAdjacents;
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

          auto firstSequence = orderedIndexSequence<3>({{
            firstAdjacents.front(),
            i,
            firstAdjacents.back()
          }});

          auto secondSequence = orderedIndexSequence<3>({{
            secondAdjacents.front(),
            i,
            secondAdjacents.back()
          }});

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
            double crossAngleLower = StdlibTypeAlgorithms::clamp(
              spiroCrossAngle(
                firstAngleBounds.upper,
                secondAngleBounds.upper
              ),
              0.0,
              M_PI
            );

            double crossAngleUpper = StdlibTypeAlgorithms::clamp(
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
                auto sequence = orderedIndexSequence<3>({{firstAdjacent, i, secondAdjacent}});

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

  /* Add default angles and dihedrals for all adjacent sequences (unpopulated
   * if e.g. the internal coordinates could not be modelled) to avoid using vdw
   * radius sums as minimum bounds for connected atom sequences.
   */
  addDefaultAngles();
  addDefaultDihedrals();
}

void SpatialModel::setBondBoundsIfEmpty(
  const std::array<AtomIndex, 2>& bondIndices,
  const double centralValue
) {
  double relativeVariance = bondRelativeVariance * _looseningMultiplier;

  /* Ensure correct use, variance for bounds should at MOST smaller than half
   * of the central value
   */
  assert(relativeVariance < 0.5 * centralValue);

  auto indexSequence = orderedIndexSequence<2>(bondIndices);

  auto findIter = _bondBounds.lower_bound(indexSequence);

  if(findIter == _bondBounds.end() || findIter->first != indexSequence) {
    _bondBounds.emplace_hint(
      findIter,
      indexSequence,
      ValueBounds {
        (1 - relativeVariance) * centralValue,
        (1 + relativeVariance) * centralValue
      }
    );
  }
}

void SpatialModel::setBondBoundsIfEmpty(
  const std::array<AtomIndex, 2>& bondIndices,
  const ValueBounds& bounds
) {
  // C++17: try_emplace
  auto indexSequence = orderedIndexSequence<2>(bondIndices);

  auto findIter = _bondBounds.lower_bound(indexSequence);

  if(findIter == _bondBounds.end() || findIter->first != indexSequence) {
    _bondBounds.emplace_hint(
      findIter,
      indexSequence,
      bounds
    );
  }
}

void SpatialModel::setAngleBoundsIfEmpty(
  const std::array<AtomIndex, 3>& angleIndices,
  const ValueBounds& bounds
) {
  auto orderedIndices = orderedIndexSequence<3>(angleIndices);

  auto findIter = _angleBounds.lower_bound(orderedIndices);

  // C++17: try_emplace
  if(findIter == _angleBounds.end() || findIter->first != orderedIndices) {
    _angleBounds.emplace_hint(
      findIter,
      orderedIndices,
      clamp(bounds, angleClampBounds)
    );
  }
}

/*!
 * Adds the dihedral bounds to the model, but only if the information for that
 * set of indices does not exist yet.
 */
void SpatialModel::setDihedralBoundsIfEmpty(
  const std::array<AtomIndex, 4>& dihedralIndices,
  const ValueBounds& bounds
) {
  auto orderedIndices = orderedIndexSequence<4>(dihedralIndices);

  auto findIter = _dihedralBounds.lower_bound(orderedIndices);

  // C++17: try_emplace
  if(findIter == _dihedralBounds.end() || findIter->first != orderedIndices) {
    _dihedralBounds.emplace_hint(
      findIter,
      orderedIndices,
      bounds
    );
  }
}

void SpatialModel::addAtomStereopermutatorInformation(
  const AtomStereopermutator& permutator,
  const std::function<double(const AtomIndex)>& cycleMultiplierForIndex,
  const double looseningMultiplier,
  const std::unordered_map<AtomIndex, Delib::Position>& fixedAngstromPositions
) {
  const auto& permutationState = permutator.getPermutationState();
  const auto& ranking = permutator.getRanking();
  const AtomIndex centerAtom = permutator.centralIndex();

  /* Angle information addition: Rough outline
   * - For each ligand
   *   - Set the distance to the center
   *     - For multi-atom ligands, correct the distance with the cone angle
   *   - Add angles between ligand-constituting atoms using the cone angle
   *     (if there are multiple ligand-constituting atoms at all)
   * - For each pair of ligands
   *   - The angle bounds between ligands are initially set by the idealized
   *     symmetry position angle modified by the upper cone angles
   *   - For each pair of atoms between both ligands
   *     - The ligand angle bounds variance is modified by cycle multipliers
   *       and loosening
   */

  /* Concerning fixed positions, if the center isn't fixed, then we do not
   * have to do anything differently at all.
   */
  std::vector<bool> ligandFixed;
  bool centerFixed = fixedAngstromPositions.count(centerAtom) > 0;
  if(centerFixed) {
    // Map from ligand indices to whether the entire ligand is fixed or not
    ligandFixed = temple::map(
      ranking.ligands,
      [&fixedAngstromPositions](const auto& ligandConstitutingIndices) -> bool {
        return temple::all_of(
          ligandConstitutingIndices,
          [&fixedAngstromPositions](const AtomIndex ligandConstitutingIndex) -> bool {
            return fixedAngstromPositions.count(ligandConstitutingIndex) > 0;
          }
        );
      }
    );
  }

  /* Intra-site modelling / For each ligand */
  for(unsigned ligandI = 0; ligandI < permutationState.ligandDistances.size(); ++ligandI) {
    /* Set the distance to the center:
     * If no cone information is present, do not correct the distance to the
     * ligand using the cone angle
     *
     * NOTE: This is probably superfluous as non-eta bonds are modelled by the
     * SpatialModel beforehand. Try removing this block below and see if that
     * causes any issues.
     */
    if(!permutationState.coneAngles.at(ligandI)) {
      for(const AtomIndex i : ranking.ligands.at(ligandI)) {
        setBondBoundsIfEmpty(
          {{i, centerAtom}},
          permutationState.ligandDistances.at(ligandI)
        );
      }

      // No further work has to be done for single-atom ligands
      continue;
    }

    if(centerFixed && ligandFixed.at(ligandI)) {
      // All center to ligand constituting atom distances are fixed *exactly*
      for(const AtomIndex i : ranking.ligands.at(ligandI)) {
        double bondDistance = DelibHelpers::distance(
          fixedAngstromPositions.at(i).asEigenVector(),
          fixedAngstromPositions.at(centerAtom).asEigenVector()
        );

        setBondBoundsIfEmpty(
          {{i, centerAtom}},
          ValueBounds {bondDistance, bondDistance}
        );
      }

      // and so are the ligand constituting atom to ligand constituting atom angles
      temple::forEach(
        temple::adaptors::allPairs(ranking.ligands.at(ligandI)),
        [&](const AtomIndex i, const AtomIndex j) {
          double angle = DelibHelpers::angle(
            fixedAngstromPositions.at(i).asEigenVector(),
            fixedAngstromPositions.at(centerAtom).asEigenVector(),
            fixedAngstromPositions.at(i).asEigenVector()
          );

          setAngleBoundsIfEmpty(
            {{i, centerAtom, j}},
            ValueBounds {angle, angle}
          );
        }
      );
    } else {
      /* Distance of every ligand site atom index to the central atom assumptions
       * - Every haptic index is on the cone base circle
       * - Cone height is defined by permutationState.ligandDistance
       * - Cone angle is defined by permutationState.coneAngle
       */
      const DistanceGeometry::ValueBounds coneAngleBounds = permutationState.coneAngles.at(ligandI).value();

      const double upperHypotenuse = (
        permutationState.ligandDistances.at(ligandI).upper
        / std::cos(coneAngleBounds.lower)
      );

      const double lowerHypotenuse = (
        permutationState.ligandDistances.at(ligandI).lower
        / std::cos(coneAngleBounds.upper)
      );

      for(const AtomIndex i : ranking.ligands.at(ligandI)) {
        setBondBoundsIfEmpty(
          {{i, centerAtom}},
          DistanceGeometry::ValueBounds {
            lowerHypotenuse,
            upperHypotenuse
          }
        );
      }

      /* Set angles between ligand-constituting atoms within a single ligand
       * - Minimally 0° (if there were a zero-length bond)
       *   You could compute shortest possible bond constexpr and insert a trig
       *   calc here, but the bond level distance is supplied elsewhere by
       *   SpatialModel anyway, no need to duplicate that information
       * - Maximally 2 * the upper cone angle (but not more than M_PI)
       */
      temple::forEach(
        temple::adaptors::allPairs(ranking.ligands.at(ligandI)),
        [&](const AtomIndex i, const AtomIndex j) {
          setAngleBoundsIfEmpty(
            {{i, centerAtom, j}},
            DistanceGeometry::ValueBounds {
              0,
              std::min(M_PI, 2 * coneAngleBounds.upper)
            }
          );
        }
      );
    }
  }

  /* Inter-site modelling / Between ligands */
  /* If for either site no cone angles could be calculated (currently only
   * happens if a haptic ligand does not match the existing modeling patterns),
   * we have to skip this step entirely and hope that the remaining modeling
   * can pick up the slack.
   *
   * NOTE: Cone angles are calculated for non-haptic ligands too -> (0, 0).
   */
  for(unsigned i = 0; i < ranking.ligands.size() - 1; ++i) {
    if(!permutationState.coneAngles.at(i)) {
      continue;
    }

    for(unsigned j = i + 1; j < ranking.ligands.size(); ++j) {
      if(!permutationState.coneAngles.at(j)) {
        continue;
      }

      if(centerFixed && ligandFixed.at(i) && ligandFixed.at(j)) {
        // All angles are known exactly!
        temple::forEach(
          temple::adaptors::allPairs(
            ranking.ligands.at(i),
            ranking.ligands.at(j)
          ),
          [&](const AtomIndex x, const AtomIndex y) -> void {
            double angle = DelibHelpers::angle(
              fixedAngstromPositions.at(x).asEigenVector(),
              fixedAngstromPositions.at(centerAtom).asEigenVector(),
              fixedAngstromPositions.at(y).asEigenVector()
            );

            setAngleBoundsIfEmpty(
              {{x, centerAtom, y}},
              ValueBounds {angle, angle}
            );
          }
        );
      } else {
        /* The idealized symmetry angles are modified by the upper (!) cone angles
         * at each site, not split between lower and upper.
         */
        ValueBounds angleBounds = makeBoundsFromCentralValue(
          permutator.angle(i, j),
          (
            permutationState.coneAngles.at(i).value().upper
            + permutationState.coneAngles.at(j).value().upper
          )
        );

        /* The computed angle bounds are valid for each pair of atoms
         * constituting each ligand
         */
        temple::forEach(
          temple::adaptors::allPairs(
            ranking.ligands.at(i),
            ranking.ligands.at(j)
          ),
          [&](const AtomIndex x, const AtomIndex y) -> void {
            double variation = (
              DistanceGeometry::SpatialModel::angleAbsoluteVariance
              * looseningMultiplier
              * cycleMultiplierForIndex(x)
              * cycleMultiplierForIndex(y)
            );

            setAngleBoundsIfEmpty(
              {{x, centerAtom, y}},
              clamp(
                ValueBounds {
                  angleBounds.lower - variation,
                  angleBounds.upper + variation
                },
                angleClampBounds
              )
            );
          }
        );
      }
    }
  }

  // Chiral constraints addition
  for(const auto& minimalConstraint : permutator.minimalChiralityConstraints()) {
    /* We need to calculate target upper and lower volumes for the chirality
     * constraints. _cache.ligandDistances contains bounds for the distance to
     * each ligand site plane, and since the center of each cone should
     * constitute the average ligand position, we can calculate 1-3 distances
     * between the centerpoints of ligands using the idealized angles.
     *
     * The target volume of the chirality constraint created by the
     * tetrahedron is calculated using internal coordinates (the
     * Cayley-Menger determinant), always leading to V > 0, so depending on
     * the current assignment, the sign of the result is switched. The
     * formula used later in chirality constraint calculation for explicit
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

    using DeterminantMatrix = Eigen::Matrix<double, 5, 5>;

    DeterminantMatrix lowerMatrix, upperMatrix;

    lowerMatrix.row(0).setOnes();
    upperMatrix.row(0).setOnes();

    lowerMatrix.diagonal().setZero();
    upperMatrix.diagonal().setZero();

    /* Cycle through all combinations of ligand indices in the tetrahedron
     * definition sequence. boost::none means the central atom.
     */
    for(unsigned i = 0; i < 4; ++i) {
      boost::optional<DistanceGeometry::ValueBounds> iBounds;
      if(minimalConstraint.at(i)) {
        iBounds = permutationState.ligandDistances.at(
          minimalConstraint.at(i).value()
        );
      }

      for(unsigned j = i + 1; j < 4; ++j) {
        boost::optional<DistanceGeometry::ValueBounds> jBounds;
        if(minimalConstraint.at(j)) {
          jBounds = permutationState.ligandDistances.at(
            minimalConstraint.at(j).value()
          );
        }

        assert(iBounds || jBounds);

        DistanceGeometry::ValueBounds oneThreeDistanceBounds;
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
              std::max(0.0, siteAngle - angleAbsoluteVariance * looseningMultiplier)
            ),
            CommonTrig::lawOfCosines(
              iBounds.value().upper,
              jBounds.value().upper,
              std::min(M_PI, siteAngle + angleAbsoluteVariance * looseningMultiplier)
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

    // Map the ligand indices to their constituent indices for use in the prototype
    auto tetrahedronLigands = temple::map(
      minimalConstraint,
      [&](const boost::optional<unsigned>& ligandIndexOptional) -> std::vector<AtomIndex> {
        if(ligandIndexOptional) {
          return ranking.ligands.at(ligandIndexOptional.value());
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
     * completely depends on the individual values. Maybe in very specific
     * cases one can deduce some relationship, but not generally.
     *
     * Helpfully, chemical_symmetry only emits positive chiral target volume
     * index sequences (see test case name allTetrahedraPositive), so no
     * negative volumes have to be considered.
     */

    _chiralConstraints.emplace_back(
      std::move(tetrahedronLigands),
      std::min(volumeFromLower, volumeFromUpper),
      std::max(volumeFromLower, volumeFromUpper)
    );
  }
}

void SpatialModel::addBondStereopermutatorInformation(
  const BondStereopermutator& permutator,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const double looseningMultiplier,
  const std::unordered_map<AtomIndex, Delib::Position>& fixedAngstromPositions
) {
  // Check preconditions and get access to commonly needed things
  assert(permutator.assigned());
  const stereopermutation::Composite& composite = permutator.composite();
  const auto& assignment = permutator.assigned();

  /* The signed volume of a tetrahedron spanned by a 1-4 sequence i-j-k-l is
   *
   *   V = 1/6 (i-l) dot [(j-l) x (k-l)]
   *
   * or, in internal coordinates (a, b, c are bond lengths, α and β angles,
   * φ the dihedral)
   *
   *   V = 1/6 a b c sin α sin β sin φ
   *
   * We calculate the maximal volume spanned by the combination of all bounds
   * on the internal coordinates and set the flat chirality constraint's
   * tolerance as ± its magnitude.
   *
   * Just like in atomStereopermutatorChiralityConstraints, we adjust the volume
   * by V' = 6 * V to avoid calculating the factor everywhere.
   *
   * As it currently stands, only dihedral approx 0° chirality constraints
   * are emitted. For more complicated composites, in which no such chirality
   * constraints exist, it would be considerably more helpful if all chirality
   * constraints were emitted.
   *
   * To change this, remove the if-condition in the loop and reconsider the
   * algorithm generating dihedral angle bounds and the maximum sine value.
   * Additionally, it will be necesssary to compute the minimum sine values
   * as well as using the lower bounds on the bounded distances to get
   * the lower volume. It will not be reasonable to clamp the dihedral to 0,
   * M_PI, but will instead be necessary to retain the +- orientation to get
   * correct +- values for the volume bounds around 0° dihedrals.
   *
   * How exactly to combine the bounds with the angle bounds and their sines to
   * get the correct bounds on the final volume will be challenging.
   */

  // Can the following selections be done with a single branch?
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

  // To determine the maximal sine value of bounds, we use boost's help
  auto sinBounds = [](const ValueBounds angleBounds) -> ValueBounds {
    using boostInterval = boost::numeric::interval<
      double,
      boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
          boost::numeric::interval_lib::rounded_transc_std<double>
        >,
        boost::numeric::interval_lib::checking_base<double>
      >
    >;

    auto transformedInterval = boost::numeric::sin(
      boostInterval(angleBounds.lower, angleBounds.upper)
    );

    return {
      boost::numeric::lower(transformedInterval),
      boost::numeric::upper(transformedInterval)
    };
  };

  /* Only dihedrals */
  unsigned firstSymmetryPosition, secondSymmetryPosition;
  double dihedralAngle;

  for(const auto& dihedralTuple : composite.dihedrals(assignment.value())) {
    std::tie(firstSymmetryPosition, secondSymmetryPosition, dihedralAngle) = dihedralTuple;

    unsigned ligandIndexIAtFirst = SymmetryMapHelper::getLigandIndexAt(
      firstSymmetryPosition,
      firstStereopermutator.getSymmetryPositionMap()
    );
    unsigned ligandIndexKAtFirst = firstStereopermutator.getRanking().getLigandIndexOf(secondStereopermutator.centralIndex());
    unsigned ligandIndexJAtSecond = secondStereopermutator.getRanking().getLigandIndexOf(firstStereopermutator.centralIndex());
    unsigned ligandIndexLAtSecond = SymmetryMapHelper::getLigandIndexAt(
      secondSymmetryPosition,
      secondStereopermutator.getSymmetryPositionMap()
    );

    const auto& coneAngleIOption = firstStereopermutator.getPermutationState().coneAngles.at(ligandIndexIAtFirst);
    const auto& coneAngleLOption = secondStereopermutator.getPermutationState().coneAngles.at(ligandIndexLAtSecond);

    // Do not emit chiral constraints if cone angles are unknown
    if(!coneAngleIOption || !coneAngleLOption) {
      continue;
    }

    const ValueBounds& coneAngleI = *coneAngleIOption;
    const ValueBounds& coneAngleL = *coneAngleLOption;

    /* Modify the dihedral angle by the upper cone angles of the i and l
     * ligands and the usual variances.
     */
    ValueBounds dihedralBounds = makeBoundsFromCentralValue(
      dihedralAngle,
      coneAngleI.upper + coneAngleL.upper + dihedralAbsoluteVariance * looseningMultiplier
    );

    /* If the width of the dihedral angle is now larger than 2π, then we may
     * overrepresent some dihedral values when choosing randomly in that
     * interval, and it is preferable to reset the interval to (-π,π].
     *
     * This should be very rare, and is just a safeguard.
     */
    if(dihedralBounds.upper - dihedralBounds.lower >= 2 * M_PI) {
      dihedralBounds = defaultDihedralBounds;
    }

    temple::forEach(
      temple::adaptors::allPairs(
        firstStereopermutator.getRanking().ligands.at(ligandIndexIAtFirst),
        secondStereopermutator.getRanking().ligands.at(ligandIndexLAtSecond)
      ),
      [&](const AtomIndex firstIndex, const AtomIndex secondIndex) -> void {
        setDihedralBoundsIfEmpty(
          std::array<AtomIndex, 4> {
            firstIndex,
            firstStereopermutator.centralIndex(),
            secondStereopermutator.centralIndex(),
            secondIndex
          },
          dihedralBounds
        );
      }
    );

    /* TODO a, c, alphaSineBounds and betaSineBounds must be modified by cone
     * angles
     */

    /* Figure out tolerance from interal bounds on bonds, angles and dihedral
     * - a is the ligand distance from the firstStereopermutator to the
     *   ligand at the symmetry position given by the composite
     * - b is the bond distance between firstStereopermutator and
     *   secondStereopermutator
     * - c is the ligand distance from the secondStereopermutator to the
     *   ligand at the symmetry position given by the composite
     * - alpha is the angle between the secondStereopermutator and the ligand
     *   at the symmetry position given by the composite at the
     *   firstStereopermutator
     * - beta is the angle between the firstStereopermutator and the ligand
     *   at the symmetry position given by the composite at the
     *   secondStereopermutator
     * - phi is the dihedralAngle (this doesn't really have any bounds yet,
     *   which we add using the dihedralAbsoluteVariance)
     */
    ValueBounds a = firstStereopermutator.getPermutationState().ligandDistances.at(ligandIndexIAtFirst);
    ValueBounds b = _bondBounds.at(
      orderedSequence(firstStereopermutator.centralIndex(), secondStereopermutator.centralIndex())
    );
    ValueBounds c = secondStereopermutator.getPermutationState().ligandDistances.at(ligandIndexLAtSecond);

    /* Construct bounds for the central values and then determine the maximal
     * sine value for each set of bounds when determining the volume
     */
    ValueBounds alphaSineBounds = sinBounds(
      clamp(
        makeBoundsFromCentralValue(
          firstStereopermutator.angle(ligandIndexIAtFirst, ligandIndexKAtFirst),
          angleAbsoluteVariance * looseningMultiplier
        ),
        angleClampBounds
      )
    );
    ValueBounds betaSineBounds = sinBounds(
      clamp(
        makeBoundsFromCentralValue(
          secondStereopermutator.angle(ligandIndexJAtSecond, ligandIndexLAtSecond),
          angleAbsoluteVariance * looseningMultiplier
        ),
        angleClampBounds
      )
    );

    /* NOTE: no clamp! Best to let the sine take care of periodicity for us
     * Also: this dihedral needs to be computed anew as it indicates the
     * dihedral between the ligand site center of masses. They are not the
     * per-atom dihedral bounds (which are modified by the ligands' cone
     * angles).
     */
    ValueBounds phiSineBounds = sinBounds(
      makeBoundsFromCentralValue(
        std::fabs(dihedralAngle),
        dihedralAbsoluteVariance * looseningMultiplier
      )
    );

    /* min and max have to be understood more in terms of most negative and
     * most positive.
     *
     * The phi sine bounds can be negative (as opposed to all other bounds)
     * and hence, if the lower phi sine bound is negative, then the most
     * negative value can be achieved with the upper value of all other bounds.
     *
     * If the lower phi sine bound is positive, then the minimal value is
     * achieved with the lower value of all bounds.
     */
    double minVolume, maxVolume;
    if(phiSineBounds.lower < 0 && phiSineBounds.upper > 0) {
      double baseVolume = a.upper * b.upper * c.upper * alphaSineBounds.upper * betaSineBounds.upper;
      minVolume = baseVolume * phiSineBounds.lower;
      maxVolume = baseVolume * phiSineBounds.upper;
    } else {
      minVolume = a.lower * b.lower * c.lower * alphaSineBounds.lower * betaSineBounds.lower * phiSineBounds.lower;
      minVolume = a.upper * b.upper * c.upper * alphaSineBounds.upper * betaSineBounds.upper * phiSineBounds.upper;
    }

    assert(minVolume < maxVolume);

    _chiralConstraints.emplace_back(
      DistanceGeometry::ChiralityConstraint::LigandSequence {
        firstStereopermutator.getRanking().ligands.at(ligandIndexIAtFirst),
        {firstStereopermutator.centralIndex()},
        {secondStereopermutator.centralIndex()},
        secondStereopermutator.getRanking().ligands.at(ligandIndexLAtSecond)
      },
      minVolume,
      maxVolume
    );
  }
}

void SpatialModel::addDefaultAngles() {
  const InnerGraph& inner = _molecule.graph().inner();
  /* If no explicit angle can be provided for a triple of bonded atoms, we need
   * to at least specify the range of possible angles so that no implicit
   * minimimum distance (sum of vdw radii) is used instead. This is important
   * for haptic ligands as for some ligand site connectivity, no circumradius
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
          {{i, center, j}},
          angleClampBounds
        );
      }
    );
  }
}

void SpatialModel::addDefaultDihedrals() {
  const InnerGraph& inner = _molecule.graph().inner();

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
            {{
              sourceAdjacentIndex,
              sourceIndex,
              targetIndex,
              targetAdjacentIndex,
            }},
            defaultDihedralBounds
          );
        }
      }
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

SpatialModel::BoundsList SpatialModel::makeBoundsList() const {
  // Copy the constraints as ground truth
  BoundsList bounds = _constraints;

  /* There may be overlapping and possibly conflicting information present in
   * the gathered data. If data affects the same atom-pair, we must ensure that
   * we merely raise the lower bound and lower the upper bound, while never
   * inverting the bounds overall.
   */
  auto addInformation = [&bounds](
    const AtomIndex i,
    const AtomIndex j,
    const ValueBounds& newBounds
  ) {
    auto indices = orderedSequence(i, j);
    auto findIter = bounds.lower_bound(indices);

    if(findIter != bounds.end() && findIter->first == indices) {
      auto& currentBounds = findIter->second;
      // Try to raise the lower bound first
      if(
        newBounds.lower > currentBounds.lower
        && newBounds.lower < currentBounds.upper
      ) {
        currentBounds.lower = newBounds.lower;
      }

      // Try to lower the upper bound
      if(
        newBounds.upper < currentBounds.upper
        && newBounds.upper > currentBounds.lower
      ) {
        currentBounds.upper = newBounds.upper;
      }
    } else {
      bounds.emplace_hint(
        findIter,
        indices,
        newBounds
      );
    }
  };

  auto getBondBounds = [&bounds](
    const AtomIndex i,
    const AtomIndex j
  ) -> ValueBounds& {
    return bounds.at(
      orderedSequence(i, j)
    );
  };

  // Add 1-2 information from the bonds
  if(_constraints.empty()) {
    /* If there are no ground constraints, then we can just copy over the
     * entire bond bounds (these are definitely compatible with triangle
     * inequalities)
     */
    bounds = _bondBounds;
  } else {
    // Otherwise, we have to carefully add bond information
    for(const auto& bondPair : _bondBounds) {
      addInformation(
        bondPair.first.front(),
        bondPair.first.back(),
        bondPair.second
      );
    }
  }

  // Add 1-3 information
  for(const auto& anglePair : _angleBounds) {
    const auto& indices = anglePair.first;
    const auto& angleBounds = anglePair.second;

    const auto& firstBounds = getBondBounds(indices.front(), indices.at(1));
    const auto& secondBounds = getBondBounds(indices.at(1), indices.back());

    addInformation(
      indices.front(),
      indices.back(),
      ValueBounds {
        CommonTrig::lawOfCosines(
          firstBounds.lower,
          secondBounds.lower,
          angleBounds.lower
        ),
        CommonTrig::lawOfCosines(
          firstBounds.upper,
          secondBounds.upper,
          angleBounds.upper
        )
      }
    );
  }

  // Add 1-4 information
  for(const auto& dihedralPair : _dihedralBounds) {
    const auto& indices = dihedralPair.first;
    const auto& dihedralBounds = dihedralPair.second;

    const auto& firstBounds = getBondBounds(indices.front(), indices.at(1));
    const auto& secondBounds = getBondBounds(indices.at(1), indices.at(2));
    const auto& thirdBounds = getBondBounds(indices.at(2), indices.back());

    auto firstAngleFindIter = _angleBounds.find(
      orderedIndexSequence<3>({{
        indices.at(0),
        indices.at(1),
        indices.at(2)
      }})
    );

    auto secondAngleFindIter = _angleBounds.find(
      orderedIndexSequence<3>({{
        indices.at(1),
        indices.at(2),
        indices.at(3)
      }})
    );

    if(
      firstAngleFindIter == _angleBounds.end()
      || secondAngleFindIter == _angleBounds.end()
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
        dihedralBounds
      )
    );
  }

  return bounds;
}

boost::optional<ValueBounds> SpatialModel::coneAngle(
  const std::vector<AtomIndex>& ligandIndices,
  const ValueBounds& coneHeightBounds
) const {
  return coneAngle(
    ligandIndices,
    coneHeightBounds,
    bondRelativeVariance * _looseningMultiplier,
    _molecule.graph(),
    Cycles {_molecule.graph(), true}
  );
}

ValueBounds SpatialModel::ligandDistance(
  const std::vector<AtomIndex>& ligandIndices,
  const AtomIndex centralIndex
) const {
  return ligandDistanceFromCenter(
    ligandIndices,
    centralIndex,
    bondRelativeVariance * _looseningMultiplier,
    _molecule.graph()
  );
}

std::vector<DistanceGeometry::ChiralityConstraint> SpatialModel::getChiralityConstraints() const {
  return _chiralConstraints;
}

void SpatialModel::dumpDebugInfo() const {
  auto& logRef = Log::log(Log::Level::Debug);
  logRef << "SpatialModel debug info" << std::endl;

  // Bonds
  for(const auto& bondIterPair : _bondBounds) {
    const auto& indexArray = bondIterPair.first;
    const auto& bounds = bondIterPair.second;
    logRef << "Bond " << temple::condense(indexArray)
      << ": [" << bounds.lower << ", " << bounds.upper << "]" << std::endl;
  }

  // Angles
  for(const auto& angleIterPair : _angleBounds) {
    const auto& indexArray = angleIterPair.first;
    const auto& bounds = angleIterPair.second;
    logRef << "Angle " << temple::condense(indexArray)
      << ": [" << bounds.lower << ", " << bounds.upper << "]" << std::endl;
  }

  // Dihedrals
  for(const auto& dihedralIterPair : _dihedralBounds) {
    const auto& indexArray = dihedralIterPair.first;
    const auto& bounds = dihedralIterPair.second;
    logRef << "Dihedral " << temple::condense(indexArray)
      << ": [" << bounds.lower << ", " << bounds.upper << "]" << std::endl;
  }
}

struct SpatialModel::ModelGraphWriter final : public MolGraphWriter {
  /* State */
  const SpatialModel& spatialModel;

/* Constructor */
  ModelGraphWriter(const InnerGraph& inner, const SpatialModel& passSpatialModel)
    : MolGraphWriter(&inner, &passSpatialModel._stereopermutators),
      spatialModel(passSpatialModel) {}

  std::vector<std::string> edgeTooltips(AtomIndex source, AtomIndex target) const final {
    const auto indexSequence = orderedIndexSequence<2>({{source, target}});
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
    std::vector<std::string> tooltips {
      Symmetry::name(permutator.getSymmetry()),
      permutator.info()
    };

    for(const auto& angleIterPair : spatialModel._angleBounds) {
      const auto& indexSequence = angleIterPair.first;
      const auto& angleBounds = angleIterPair.second;

      if(indexSequence.at(1) == permutator.centralIndex()) {
        tooltips.emplace_back(
          "["s + std::to_string(indexSequence.at(0)) + ","s
          + std::to_string(indexSequence.at(2)) +"] -> ["s
          + std::to_string(
            temple::Math::round(
              temple::Math::toDegrees(angleBounds.lower)
            )
          ) + ", "s + std::to_string(
            temple::Math::round(
              temple::Math::toDegrees(angleBounds.upper)
            )
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
            temple::Math::round(
              temple::Math::toDegrees(dihedralBounds.lower)
            )
          ) + ", "s + std::to_string(
            temple::Math::round(
              temple::Math::toDegrees(dihedralBounds.upper)
            )
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
boost::optional<ValueBounds> SpatialModel::coneAngle(
  const std::vector<AtomIndex>& baseConstituents,
  const ValueBounds& coneHeightBounds,
  const double bondRelativeVariance,
  const OuterGraph& graph,
  const Cycles& etaLessCycles
) {
  /* Have to decide cone base radius in order to calculate this. There are some
   * simple cases to get out of the way first:
   */

  assert(!baseConstituents.empty());
  if(baseConstituents.size() == 1) {
    return ValueBounds {0.0, 0.0};
  }

  if(baseConstituents.size() == 2) {
    double radius = Bond::calculateBondDistance(
      graph.elementType(baseConstituents.front()),
      graph.elementType(baseConstituents.back()),
      graph.bondType(
        BondIndex {
          baseConstituents.front(),
          baseConstituents.back()
        }
      )
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
  /* This is essentially an if - only one cycle can consist of exactly the base
   * constituents
   */
  for(
    const auto cyclePtr :
    boost::make_iterator_range(
      etaLessCycles.iteratorPair(Cycles::predicates::ConsistsOf {baseConstituents})
    )
  ) {
    /* So if it IS a cycle, we need a ring index sequence to calculate a cyclic
     * polygon circumradius, which is how flat cycles are modelled here
     */
    auto ringIndexSequence = makeRingIndexSequence(
      Cycles::edgeVertices(cyclePtr)
    );

    auto distances = temple::map(
      temple::adaptors::sequentialPairs(ringIndexSequence),
      [&](const AtomIndex i, const AtomIndex j) -> double {
        return Bond::calculateBondDistance(
          graph.elementType(i),
          graph.elementType(j),
          graph.bondType(BondIndex {i, j})
        );
      }
    );

    auto lowerCircumradiusResult = CyclicPolygons::detail::convexCircumradius(
      temple::map(
        distances,
        [&](const double distance) -> double {
          return (1 - bondRelativeVariance) * distance;
        }
      )
    );

    auto upperCircumradiusResult = CyclicPolygons::detail::convexCircumradius(
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


  /* So the ligand atoms are NOT the sole constituents of a closed cycle.
   *
   * For some types of ligands, we could still figure out a cone angle. If the
   * ligand group is actually a path in which any intermediate atom merely
   * connects the subsequent atoms (i.e. the path is not branched, there are no
   * cycles), and if all the involved intermediate geometries constituting the
   * longest path have only one distinct angle value, then we can create a
   * conformational model anyway.
   *
   * However, a path specific approach cannot treat branched haptic ligands
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

ValueBounds SpatialModel::ligandDistanceFromCenter(
  const std::vector<AtomIndex>& ligandIndices,
  const AtomIndex centralIndex,
  const double bondRelativeVariance,
  const OuterGraph& graph
) {
  assert(!ligandIndices.empty());

  Delib::ElementType centralIndexType = graph.elementType(centralIndex);

  if(ligandIndices.size() == 1) {
    auto ligandIndex = ligandIndices.front();

    double distance = Bond::calculateBondDistance(
      graph.elementType(ligandIndex),
      centralIndexType,
      graph.bondType(
        BondIndex {ligandIndex, centralIndex}
      )
    );

    return {
      (1 - bondRelativeVariance) * distance,
      (1 + bondRelativeVariance) * distance
    };
  }

  double distance = 0.9 * temple::average(
    temple::adaptors::transform(
      ligandIndices,
      [&](AtomIndex ligandIndex) -> double {
        return Bond::calculateBondDistance(
          graph.elementType(ligandIndex),
          centralIndexType,
          graph.bondType(
            BondIndex {ligandIndex, centralIndex}
          )
        );
      }
    )
  );

  return {
    (1 - bondRelativeVariance) * distance,
    (1 + bondRelativeVariance) * distance
  };
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
  const ValueBounds& bounds,
  const ValueBounds& clampBounds
) {
  return {
    StdlibTypeAlgorithms::clamp(
      bounds.lower,
      clampBounds.lower,
      clampBounds.upper
    ),
    StdlibTypeAlgorithms::clamp(
      bounds.upper,
      clampBounds.lower,
      clampBounds.upper
    )
  };
}

void SpatialModel::checkFixedPositionsPreconditions(
  const Molecule& molecule,
  const Configuration& configuration
) {
  // Early exit if there are no fixed positions
  if(configuration.fixedPositions.empty()) {
    return;
  }

  /* Check to ensure that every fixed atom has zero, one or all ligand sites
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
      // Check to ensure either 0, 1 or L ligand sites are fixed
      unsigned numFixedLigands = temple::accumulate(
        stereopermutatorOption->getRanking().ligands,
        0u,
        [&fixedAtoms](const unsigned carry, const auto& indexSet) -> unsigned {
          const unsigned countFixed = temple::accumulate(
            indexSet,
            0u,
            [&fixedAtoms](const unsigned nestedCarry, const AtomIndex i) -> unsigned {
              return nestedCarry + fixedAtoms.count(i);
            }
          );

          if(countFixed != 0 && countFixed != indexSet.size()) {
            throw std::runtime_error(
              "DG preconditions for fixed atoms are not met: A non-terminal "
              "atom's ligand site is only partially fixed."
            );
          }

          // Count this ligand if the site is fully fixed
          if(countFixed == indexSet.size()) {
            return carry + 1;
          }

          return carry;
        }
      );

      if(1 < numFixedLigands && numFixedLigands < Symmetry::size(stereopermutatorOption->getSymmetry())) {
        throw std::runtime_error(
          "DG preconditions for fixed atoms are not met: A non-terminal atom "
          "does not have 0, 1 or all ligand sites fixed."
        );
      }
    }
  }
}

} // namespace DistanceGeometry

} // namespace molassembler
