#include "molassembler/DistanceGeometry/SpatialModel.h"

#include "boost/graph/graphviz.hpp"
#include "Delib/ElementInfo.h"
#include "CyclicPolygons.h"

#include "chemical_symmetries/Properties.h"
#include "temple/Random.h"
#include "temple/Containers.h"

#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/detail/CommonTrig.h"
#include "molassembler/detail/StdlibTypeAlgorithms.h"
#include "molassembler/detail/MolGraphWriter.h"
#include "molassembler/Cycles.h"
#include "molassembler/Log.h"
#include "molassembler/RankingInformation.h"

#include <fstream>

namespace molassembler {

namespace DistanceGeometry {

// General availability of static constexpr members
constexpr double SpatialModel::bondRelativeVariance;
constexpr double SpatialModel::angleAbsoluteVariance;
constexpr double SpatialModel::dihedralAbsoluteVariance;

const ValueBounds SpatialModel::angleClampBounds = ValueBounds {0.0, M_PI};
const ValueBounds SpatialModel::dihedralClampBounds = ValueBounds {0.0, M_PI};

SpatialModel::SpatialModel(
  const Molecule& molecule,
  const double looseningMultiplier
) : _molecule(molecule),
    _looseningMultiplier(looseningMultiplier),
    _stereocenters(molecule.getStereocenterList())
{
  /* This is overall a pretty complicated constructor since it encompasses the
   * entire conversion from a molecular graph into some model of the relations
   * between atoms, determining which conformations are accessible.
   *
   * The rough sequence of operations is
   * - Make helper variables
   * - Set 1-2 bounds.
   * - Gather information on local geometries of all non-terminal atoms, using
   *   the existing stereocenter data and supplanting it with random-assignment
   *   inferred stereocenters on all other non-terminal atoms or double bonds.
   *
   *   - Copy original molecule's list of stereocenters. Set unassigned ones to
   *     a random assignment (consistent with occurrence statistics)
   *   - Instantiate BondStereocenters on all double bonds that aren't immediately
   *     involved in other stereocenters. This is to ensure that all atoms
   *     involved in non-stereogenic double bonds are also flat in the
   *     resulting 3D structure. They additionally allow extraction of
   *     angle and dihedral angle information just like AtomStereocenters.
   *   - Instantiate AtomStereocenters on all remaining atoms and default-assign
   *     them. This is so that we can get angle data between substituents easily.
   *
   * - Set internal angles of all small flat cycles
   * - Set all remaining 1-3 bounds with additional tolerance if atoms involved
   *   in the angle are part of a small cycle
   * - Add BondStereocenter 1-4 bound information
   */

  // Helper variables
  // Ignore eta bonds in construction of cycle data
  Cycles cycleData {molecule.getGraph(), true};
  //Cycles cycleData = molecule.getCycleData();
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

  // Set bond distances
  for(const auto& edge: molecule.iterateEdges()) {
    auto bondType = molecule.getGraph()[edge].bondType;

    // Do not model eta bonds, stereocenters are responsible for those
    if(bondType == BondType::Eta) {
      continue;
    }

    AtomIndexType i = boost::source(edge, molecule.getGraph());
    AtomIndexType j = boost::target(edge, molecule.getGraph());

    auto bondDistance = Bond::calculateBondDistance(
      molecule.getElementType(i),
      molecule.getElementType(j),
      bondType
    );

    setBondBoundsIfEmpty(
      {{i, j}},
      bondDistance
    );
  }

  /* The StereocenterList is already copy initialized with the Molecule's
   * stereocenters, but we need to instantiate AtomStereocenters everywhere,
   * regardless of whether they are stereogenic or not, to ensure that
   * modelling gets the information it needs.
   *
   * So for every missing non-terminal atom, create a AtomStereocenter in the
   * determined geometry
   */
  for(unsigned i = 0; i < molecule.numAtoms(); ++i) {
    // Already an instantiated AtomStereocenter?
    if(_stereocenters.option(i)) {
      continue;
    }

    auto localRanking = molecule.rankPriority(i);

    // Terminal atom index?
    if(localRanking.ligands.size() <= 1) {
      continue;
    }

    Symmetry::Name localSymmetry = molecule.determineLocalGeometry(i, localRanking);

    auto newStereocenter = AtomStereocenter {
      molecule.getGraph(),
      localSymmetry,
      i,
      std::move(localRanking)
    };

    /* New stereocenters encountered at this point can have multiple
     * assignments, since some types of stereocenters are flatly ignored by the
     * candidate functions from Molecule, such as trigonal pyramidal nitrogens.
     * These are found here, though, and MUST be chosen randomly according to
     * the relative weights to get a single conformation in the final model
     */
    newStereocenter.assignRandom();

    // Add it to the list of stereocenters
    _stereocenters.add(
      i,
      std::move(newStereocenter)
    );
  }

  /* For all flat cycles for which the internal angles can be determined
   * exactly, add that information
   */
  for(
    const auto cyclePtr :
    cycleData.iterate(Cycles::predicates::SizeLessThan {6})
  ) {
    const auto edgeDescriptors = Cycles::edges(cyclePtr, molecule.getGraph());
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
          molecule.getGraph()
        ) >= 1
      )
      /* TODO missing cases:
       * - Need aromaticity checking routine for cycle size 5
       */
    ) {
      /* Gather sequence of atoms in cycle by progressively converting edge
       * descriptors into vertex indices
       */
      const auto indexSequence = makeRingIndexSequence(
        Cycles::edgeVertices(cyclePtr)
      );

      /* First, we fetch the angles that maximize the cycle area using the
       * cyclic polygon library.
       */
      const auto cycleInternalAngles = CyclicPolygons::internalAngles(
        // Map sequential index pairs to their purported bond length
        temple::mapSequentialPairs(
          indexSequence,
          [&](const AtomIndexType i, const AtomIndexType j) -> double {
            auto bondTypeOption = molecule.getBondType(i, j);

            // These vertices really ought to be bonded
            assert(bondTypeOption);

            return Bond::calculateBondDistance(
              molecule.getElementType(i),
              molecule.getElementType(j),
              bondTypeOption.value()
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
  auto cycleMultiplierForIndex = [&](const AtomIndexType i) -> double {
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

  // Get 1-3 information from AtomStereocenters
  for(const auto& stereocenter : _stereocenters.atomStereocenters()) {
    stereocenter.setModelInformation(
      *this,
      cycleMultiplierForIndex,
      looseningMultiplier
    );
  }

  // Get 1-4 information from BondStereocenters
  for(const auto& bondStereocenter : _stereocenters.bondStereocenters()) {
    bondStereocenter.setModelInformation(
      *this,
      _stereocenters.option(
        boost::source(bondStereocenter.edge(), molecule.getGraph())
      ).value(),
      _stereocenters.option(
        boost::target(bondStereocenter.edge(), molecule.getGraph())
      ).value(),
      looseningMultiplier
    );
  }

  /* Model spiro centers */
  /* For an atom to be a spiro center, it needs to be contained in exactly two
   * URFs and have four substituents in a tetrahedral symmetry.
   *
   * We also only have to worry about modeling them if the two URFs are
   * particularly small, i.e. are sizes 3-5
   */
  for(const auto& stereocenter : _stereocenters.atomStereocenters()) {
    AtomIndexType i = stereocenter.centralIndex();

    // Skip any stereocenters that do not match our conditions
    if(
      stereocenter.getSymmetry() != Symmetry::Name::Tetrahedral
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
        auto makeVerticesSet = [](const auto& cyclePtr) -> std::set<AtomIndexType> {
          std::set<AtomIndexType> vertices;

          for(unsigned i = 0; i < cyclePtr -> weight; ++i) {
            vertices.insert(cyclePtr->edges[i][0]);
            vertices.insert(cyclePtr->edges[i][1]);
          }

          return vertices;
        };

        auto cycleOneVertices = makeVerticesSet(cycleOne);
        auto cycleTwoVertices = makeVerticesSet(cycleTwo);

        auto intersection = temple::setIntersection(
          cycleOneVertices,
          cycleTwoVertices
        );

        if(intersection.size() == 1 && *intersection.begin() == i) {
          auto iAdjacents = molecule.getAdjacencies(i);

          std::vector<AtomIndexType> firstAdjacents, secondAdjacents;
          for(const auto& iAdjacent : iAdjacents) {
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

            temple::forAllPairs(
              firstAdjacents,
              secondAdjacents,
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
                  // TODO maybe log this behavior?
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

  // Add default angles and dihedrals for all adjacent sequences
  addDefaultAngles();
  addDefaultDihedrals();
}

void SpatialModel::setBondBoundsIfEmpty(
  const std::array<AtomIndexType, 2>& bondIndices,
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
  const std::array<AtomIndexType, 2>& bondIndices,
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
  const std::array<AtomIndexType, 3>& angleIndices,
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
  const std::array<AtomIndexType, 4>& dihedralIndices,
  const ValueBounds& bounds
) {
  auto orderedIndices = orderedIndexSequence<4>(dihedralIndices);

  auto findIter = _dihedralBounds.lower_bound(orderedIndices);

  // C++17: try_emplace
  if(findIter == _dihedralBounds.end() || findIter->first != orderedIndices) {
    _dihedralBounds.emplace_hint(
      findIter,
      orderedIndices,
      clamp(bounds, dihedralClampBounds)
    );
  }
}

void SpatialModel::addDefaultAngles() {
  /* If no explicit angle can be provided for a triple of bonded atoms, we need
   * to at least specify the range of possible angles so that no implicit
   * minimimum distance (sum of vdw radii) is used instead.
   */

  const AtomIndexType N = _molecule.numAtoms();
  for(AtomIndexType center = 0; center < N; ++center) {
    temple::forAllPairs(
      _molecule.getAdjacencies(center),
      [&](const AtomIndexType i, const AtomIndexType j) -> void {
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
  for(const auto& edgeDescriptor : _molecule.iterateEdges()) {
    const AtomIndexType sourceIndex = boost::source(edgeDescriptor, _molecule.getGraph());
    const AtomIndexType targetIndex = boost::target(edgeDescriptor, _molecule.getGraph());

    auto sourceAdjacencies = _molecule.getAdjacencies(sourceIndex);
    auto targetAdjacencies = _molecule.getAdjacencies(targetIndex);

    temple::inplaceRemove(sourceAdjacencies, targetIndex);
    temple::inplaceRemove(targetAdjacencies, sourceIndex);

    temple::forAllPairs(
      sourceAdjacencies,
      targetAdjacencies,
      [&](
        const AtomIndexType sourceAdjacentIndex,
        const AtomIndexType targetAdjacentIndex
      ) -> void {
        if(sourceAdjacentIndex != targetAdjacentIndex) {
          setDihedralBoundsIfEmpty(
            {{
              sourceAdjacentIndex,
              sourceIndex,
              targetIndex,
              targetAdjacentIndex,
            }},
            dihedralClampBounds
          );
        }
      }
    );
  }
}

template<std::size_t N>
bool bondInformationIsPresent(
  const DistanceBoundsMatrix& bounds,
  const std::array<AtomIndexType, N>& indices
) {
  // Ensure all indices are unique
  std::set<AtomIndexType> indicesSet {indices.begin(), indices.end()};
  if(indicesSet.size() < indices.size()) {
    return false;
  }

  // Check that the bond information in the sequence is present
  return temple::all_of(
    temple::mapSequentialPairs(
      indices,
      [&bounds](const AtomIndexType i, const AtomIndexType j) -> bool {
        return (
          bounds.lowerBound(i, j) != DistanceBoundsMatrix::defaultLower
          && bounds.upperBound(i, j) != DistanceBoundsMatrix::defaultUpper
        );
      }
    )
  );
}

SpatialModel::BoundsList SpatialModel::makeBoundsList() const {
  BoundsList bounds = _bondBounds;

  auto addInformation = [&bounds](
    const AtomIndexType i,
    const AtomIndexType j,
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
    const AtomIndexType i,
    const AtomIndexType j
  ) -> ValueBounds& {
    return bounds.at(
      orderedSequence(i, j)
    );
  };

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
      ValueBounds {
        CommonTrig::dihedralLength(
          firstBounds.lower,
          secondBounds.lower,
          thirdBounds.lower,
          abAngleBounds.lower,
          bcAngleBounds.lower,
          dihedralBounds.lower // cis dihedral
        ),
        CommonTrig::dihedralLength(
          firstBounds.upper,
          secondBounds.upper,
          thirdBounds.upper,
          abAngleBounds.upper,
          bcAngleBounds.upper,
          dihedralBounds.upper // trans dihedral
        )
      }
    );
  }

  return bounds;
}

boost::optional<ValueBounds> SpatialModel::coneAngle(
  const std::vector<AtomIndexType>& ligandIndices,
  const ValueBounds& coneHeightBounds
) const {
  return coneAngle(
    ligandIndices,
    coneHeightBounds,
    bondRelativeVariance * _looseningMultiplier,
    _molecule.getGraph(),
    Cycles {_molecule.getGraph(), true}
  );
}

ValueBounds SpatialModel::ligandDistance(
  const std::vector<AtomIndexType>& ligandIndices,
  const AtomIndexType centralIndex
) const {
  return ligandDistanceFromCenter(
    ligandIndices,
    centralIndex,
    bondRelativeVariance * _looseningMultiplier,
    _molecule.getGraph()
  );
}

std::vector<DistanceGeometry::ChiralityConstraint> SpatialModel::getChiralityConstraints() const {
  std::vector<DistanceGeometry::ChiralityConstraint> constraints;

  for(const auto& stereocenter : _stereocenters.atomStereocenters()) {
    auto constraints = stereocenter.chiralityConstraints(_looseningMultiplier);

    std::move(
      std::begin(constraints),
      std::end(constraints),
      std::back_inserter(constraints)
    );
  }

  for(const auto& bondStereocenter : _stereocenters.bondStereocenters()) {
    auto constraints = bondStereocenter.chiralityConstraints(
      _looseningMultiplier,
      _stereocenters.option(
        boost::source(bondStereocenter.edge(), _molecule.getGraph())
      ).value(),
      _stereocenters.option(
        boost::target(bondStereocenter.edge(), _molecule.getGraph())
      ).value()
    );

    std::move(
      std::begin(constraints),
      std::end(constraints),
      std::back_inserter(constraints)
    );
  }

  return constraints;
}

void SpatialModel::dumpDebugInfo() const {
  auto& logRef = Log::log(Log::Level::Debug);
  logRef << "SpatialModel debug info" << std::endl;

  // Bonds
  for(const auto& bondIterPair : _bondBounds) {
    const auto& indexArray = bondIterPair.first;
    const auto& bounds = bondIterPair.second;
    logRef << "Bond " << temple::condenseIterable(indexArray)
      << ": [" << bounds.lower << ", " << bounds.upper << "]" << std::endl;
  }

  // Angles
  for(const auto& angleIterPair : _angleBounds) {
    const auto& indexArray = angleIterPair.first;
    const auto& bounds = angleIterPair.second;
    logRef << "Angle " << temple::condenseIterable(indexArray)
      << ": [" << bounds.lower << ", " << bounds.upper << "]" << std::endl;
  }

  // Dihedrals
  for(const auto& dihedralIterPair : _dihedralBounds) {
    const auto& indexArray = dihedralIterPair.first;
    const auto& bounds = dihedralIterPair.second;
    logRef << "Dihedral " << temple::condenseIterable(indexArray)
      << ": [" << bounds.lower << ", " << bounds.upper << "]" << std::endl;
  }
}

struct SpatialModel::ModelGraphWriter {
  /* State */
  // We promise to be good and not change anything
  const GraphType* const graphPtr;
  const SpatialModel& spatialModel;

/* Constructor */
  ModelGraphWriter(
    const GraphType* passGraphPtr,
    const SpatialModel& spatialModel
  ) : graphPtr(passGraphPtr),
    spatialModel(spatialModel)
  {}

/* Helper functions */
  Delib::ElementType getElementType(const AtomIndexType& vertexIndex) const {
    return (*graphPtr)[vertexIndex].elementType;
  }

/* Accessors for boost::write_graph */
  // Global options
  void operator() (std::ostream& os) const {
    os << "  graph [fontname = \"Arial\", layout = neato];\n"
      << "  node [fontname = \"Arial\", shape = circle, style = filled];\n"
      << "  edge [fontname = \"Arial\"];\n";

    /* Additional nodes */
    for(const auto& bondStereocenter : spatialModel._stereocenters.bondStereocenters()) {
      std::string state;
      if(bondStereocenter.assigned()) {
        state = std::to_string(bondStereocenter.assigned().value());
      } else {
        state = "u";
      }

      state += "/"s + std::to_string(bondStereocenter.numStereopermutations());

      std::string graphNodeName = "BS" + temple::condenseIterable(
        spatialModel._molecule.vertices(bondStereocenter.edge()),
        ""
      );

      os << "  " << graphNodeName << R"( [label=")" << state
        << R"(", fillcolor="steelblue", shape="box", fontcolor="white", )"
        << R"(tooltip=")"
        << bondStereocenter.info()
        << R"("];)" << "\n";

      // Add connections to the vertices (although those don't exist yet)
      for(const auto& vertex : spatialModel._molecule.vertices(bondStereocenter.edge())) {
        os << "  " << graphNodeName << " -- " << vertex
          << R"( [color="gray", dir="forward", len="2"];)"
          << "\n";
      }
    }
  }

  // Vertex options
  void operator() (std::ostream& os, const AtomIndexType& vertexIndex) const {
    const std::string symbolString = Delib::ElementInfo::symbol(
      getElementType(vertexIndex)
    );

    os << "[";

    // Add element name and index label
    os << R"(label = ")" << symbolString << vertexIndex << R"(")";

    // Coloring
    // C++17 if-init
    auto bgColorFindIter = MolGraphWriter::elementBGColorMap.find(symbolString);
    if(bgColorFindIter != MolGraphWriter::elementBGColorMap.end()) {
      os << R"(, fillcolor=")" << bgColorFindIter->second << R"(")";
    } else { // default
      os << R"(, fillcolor="white")";
    }

    auto textColorFindIter = MolGraphWriter::elementTextColorMap.find(symbolString);
    if(textColorFindIter != MolGraphWriter::elementTextColorMap.end()) {
      os << R"(, fontcolor=")" << textColorFindIter->second << R"(")";
    } else { // default
      os << R"(, fontcolor="orange")";
    }

    // Font sizing
    if(symbolString == "H") {
      os << ", fontsize=10, width=.3, fixedsize=true";
    }

    // Any angles this atom is the central atom in
    std::vector<std::string> tooltipStrings;

    if(auto stereocenterOption = spatialModel._stereocenters.option(vertexIndex)) {
      tooltipStrings.emplace_back(
        Symmetry::name(stereocenterOption.value().getSymmetry())
      );
      tooltipStrings.emplace_back(
        stereocenterOption.value().info()
      );
    }

    for(const auto& angleIterPair : spatialModel._angleBounds) {
      const auto& indexSequence = angleIterPair.first;
      const auto& angleBounds = angleIterPair.second;

      if(indexSequence.at(1) == vertexIndex) {
        tooltipStrings.emplace_back(
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

    if(!tooltipStrings.empty()) {
      os << R"(, tooltip=")"
        << temple::condenseIterable(tooltipStrings, "&#10;"s)
        << R"(")";
    }

    os << "]";
  }

  // Edge options
  void operator() (std::ostream& os, const EdgeIndexType& edgeIndex) const {
    const AtomIndexType source = boost::source(edgeIndex, *graphPtr);
    const AtomIndexType target = boost::target(edgeIndex, *graphPtr);

    os << "[";

    // Bond Type display options
    auto bondType = (*graphPtr)[edgeIndex].bondType;
    auto stringFindIter = MolGraphWriter::bondTypeDisplayString.find(bondType);
    if(stringFindIter != MolGraphWriter::bondTypeDisplayString.end()) {
      os << stringFindIter->second;
    }

    // If one of the bonded atoms is a hydrogen, shorten the bond
    /*if(
      getElementType(target) == Delib::ElementType::H
      || getElementType(source) == Delib::ElementType::H
    ) {
      os << ", len=0.5";
    }*/

    os << ", penwidth=3";

    const auto indexSequence = orderedIndexSequence<2>({{source, target}});
    if(spatialModel._bondBounds.count(indexSequence) == 1) {
      const auto& bondBounds = spatialModel._bondBounds.at(indexSequence);
      os << R"(, edgetooltip="[)" << bondBounds.lower << ", " << bondBounds.upper <<  R"(]")";
    }

    os << "]";
  }
};

std::string SpatialModel::dumpGraphviz() const {
  ModelGraphWriter graphWriter(
    &_molecule.getGraph(),
    *this
  );

  std::stringstream graphvizStream;

  boost::write_graphviz(
    graphvizStream,
    _molecule.getGraph(),
    graphWriter,
    graphWriter,
    graphWriter
  );

  return graphvizStream.str();
}

void SpatialModel::writeGraphviz(const std::string& filename) const {
  ModelGraphWriter graphWriter(
    &_molecule.getGraph(),
    *this
  );

  std::ofstream outStream(filename);

  boost::write_graphviz(
    outStream,
    _molecule.getGraph(),
    graphWriter,
    graphWriter,
    graphWriter
  );

  outStream.close();
}

/* Static functions */
boost::optional<ValueBounds> SpatialModel::coneAngle(
  const std::vector<AtomIndexType>& baseConstituents,
  const ValueBounds& coneHeightBounds,
  const double bondRelativeVariance,
  const GraphType& graph,
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
    auto findEdgePair = boost::edge(
      baseConstituents.front(),
      baseConstituents.back(),
      graph
    );

    // If a ligand is haptic and consists of two ligands, these MUST be bonded
    assert(findEdgePair.second);

    double radius = Bond::calculateBondDistance(
      graph[baseConstituents.front()].elementType,
      graph[baseConstituents.back()].elementType,
      graph[findEdgePair.first].bondType
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
    etaLessCycles.iterate(Cycles::predicates::ConsistsOf {baseConstituents})
  ) {
    /* So if it IS a cycle, we need a ring index sequence to calculate a cyclic
     * polygon circumradius, which is how flat cycles are modelled here
     */
    auto ringIndexSequence = makeRingIndexSequence(
      Cycles::edgeVertices(cyclePtr)
    );

    auto distances = temple::mapSequentialPairs(
      ringIndexSequence,
      [&](const AtomIndexType i, const AtomIndexType j) -> double {
        auto findEdgePair = boost::edge(i, j, graph);

        // These vertices really ought to be bonded
        assert(findEdgePair.second);

        return Bond::calculateBondDistance(
          graph[i].elementType,
          graph[j].elementType,
          graph[findEdgePair.first].bondType
        );
      }
    );

    auto lowerCircumradiusResult = CyclicPolygons::detail::convexCircumradius(
      temple::map(
        distances,
        [&](const double& distance) -> double {
          return (1 - bondRelativeVariance) * distance;
        }
      )
    );

    auto upperCircumradiusResult = CyclicPolygons::detail::convexCircumradius(
      temple::map(
        distances,
        [&](const double& distance) -> double {
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
   * to the Molecule's StereocenterList in both cases. In that case, this
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
  const std::vector<AtomIndexType>& ligandIndices,
  const AtomIndexType centralIndex,
  const double bondRelativeVariance,
  const GraphType& graph
) {
  assert(!ligandIndices.empty());

  Delib::ElementType centralIndexType = graph[centralIndex].elementType;

  if(ligandIndices.size() == 1) {
    auto ligandIndex = ligandIndices.front();

    auto edgeFindPair = boost::edge(ligandIndex, centralIndex, graph);
    assert(edgeFindPair.second);

    double distance = Bond::calculateBondDistance(
      graph[ligandIndex].elementType,
      centralIndexType,
      graph[edgeFindPair.first].bondType
    );

    return {
      (1 - bondRelativeVariance) * distance,
      (1 + bondRelativeVariance) * distance
    };
  }

  double distance = 0.9 * temple::average(
    temple::map(
      ligandIndices,
      [&](AtomIndexType ligandIndex) -> double {
        auto edgeFindPair = boost::edge(ligandIndex, centralIndex, graph);
        assert(edgeFindPair.second);

        return Bond::calculateBondDistance(
          graph[ligandIndex].elementType,
          centralIndexType,
          graph[edgeFindPair.first].bondType
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

} // namespace DistanceGeometry

} // namespace molassembler
