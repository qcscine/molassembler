#include "boost/graph/graphviz.hpp"
#include "CyclicPolygons.h"
#include "Delib/ElementInfo.h"
#include "chemical_symmetries/Properties.h"

#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "DistanceGeometry/DistanceGeometry.h"
#include "CommonTrig.h"
#include "Log.h"
#include "StdlibTypeAlgorithms.h"
#include "temple/Random.h"
#include "temple/Containers.h"

#include <fstream>

/* TODO
 * - Use the static const MolGraphWriter coloring maps
 */

namespace molassembler {

namespace DistanceGeometry {

// General availability of static constexpr members
constexpr double MoleculeSpatialModel::bondRelativeVariance;
constexpr double MoleculeSpatialModel::angleAbsoluteVariance;
constexpr double MoleculeSpatialModel::dihedralAbsoluteVariance;

MoleculeSpatialModel::MoleculeSpatialModel(
  const Molecule& molecule,
  const double looseningMultiplier
) : _molecule(molecule), _looseningMultiplier(looseningMultiplier) {
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
   *   - Instantiate EZStereocenters on all double bonds that aren't immediately
   *     involved in other stereocenters. This is to ensure that all atoms
   *     involved in non-stereogenic double bonds are also flat in the
   *     resulting 3D structure. They additionally allow extraction of
   *     angle and dihedral angle information just like CNStereocenters.
   *   - Instantiate CNStereocenters on all remaining atoms and default-assign
   *     them. This is so that we can get angle data between substituents easily.
   *
   * - Set internal angles of all small flat cycles
   * - Set all remaining 1-3 bounds with additional tolerance if atoms involved
   *   in the angle are part of a small cycle
   * - Add EZStereocenter 1-4 bound information
   */

  // Helper variables
  // Ignore eta bonds in construction of cycle data
  Cycles cycleData {molecule.getGraph(), true};
  //Cycles cycleData = molecule.getCycleData();
  auto smallestCycleMap = makeSmallestCycleMap(cycleData);

  // Check constraints on static constants
  static_assert(
    0 < bondRelativeVariance && bondRelativeVariance < 1,
    "MoleculeSpatialModel static constant bond relative variance must fulfill"
    "0 < x << 1"
  );
  static_assert(
    0 < angleAbsoluteVariance && angleAbsoluteVariance < Symmetry::smallestAngle,
    "MoleculeSpatialModel static constant angle absolute variance must fulfill"
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

  // Populate the stereocenterMap with copies of the molecule's stereocenters
  // Start with the existing stereocenters with multiple assignments
  for(const auto& stereocenterPtr : molecule.getStereocenterList()) {
    for(const auto& involvedAtom : stereocenterPtr -> involvedAtoms()) {
      if(_stereocenterMap.count(involvedAtom) == 0) {
        if(stereocenterPtr -> type() == Stereocenters::Type::CNStereocenter) {
          // Downcast the shared ptr
          auto CNSPtr = std::dynamic_pointer_cast<
            Stereocenters::CNStereocenter
          >(stereocenterPtr);

          // Explicit new shared ptr using derived copy-constructor
          _stereocenterMap[involvedAtom] = std::make_shared<
            Stereocenters::CNStereocenter
          >(*CNSPtr);

          /* Now we have a stereocenter, but it might be unassigned, in which
           * case angle() will fail! Need to assign unassigned ones
           * at random consistent with the unique assignments' relative
           * occurrences
           */
          if(!_stereocenterMap[involvedAtom] -> assigned()) {
            std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
              _stereocenterMap[involvedAtom]
            ) -> assignRandom();
          }
        } else {
          auto EZSPtr = std::dynamic_pointer_cast<
            Stereocenters::EZStereocenter
          >(stereocenterPtr);

          _stereocenterMap[involvedAtom] = std::make_shared<
            Stereocenters::EZStereocenter
          >(*EZSPtr);

          if(!_stereocenterMap[involvedAtom] -> assigned()) {
            // Assign the EZStereocenter at random
            _stereocenterMap[involvedAtom] -> assign(
              temple::random.getSingle<unsigned>(
                0,
                _stereocenterMap[involvedAtom] -> numStereopermutations() - 1
              )
            );
          }
        }
      }
    }
  }

  /* For every non-implicated double bond, create an EZStereocenter
   * Implicated means here that no stereocenters exist for either of the double
   * edge's atoms.
   *
   * TODO awkwardness here: double / aromatic bond type schism, aromatic bond
   * type not considered
   */
  for(const auto& edgeIndex: molecule.iterateEdges()) {
    if(molecule.getGraph()[edgeIndex].bondType == BondType::Double) {
      auto source = boost::source(edgeIndex, molecule.getGraph()),
           target = boost::target(edgeIndex, molecule.getGraph());

      if(
        _stereocenterMap.count(source) == 0
        && _stereocenterMap.count(target) == 0
        && molecule.getNumAdjacencies(source) == 3
        && molecule.getNumAdjacencies(target) == 3
      ) {
        // Instantiate without regard for number of assignments
        auto newStereocenterPtr = std::make_shared<Stereocenters::EZStereocenter>(
          source,
          molecule.rankPriority(source, {target}), // exclude shared edge
          target,
          molecule.rankPriority(target, {source})
        );

        // Map source *and* target to the same stereocenterPtr
        _stereocenterMap[source] = newStereocenterPtr;
        _stereocenterMap[target] = newStereocenterPtr;
      }
    }
  }

  /* For every missing non-terminal atom, create a CNStereocenter in the
   * determined geometry
   */
  for(unsigned i = 0; i < molecule.numAtoms(); i++) {
    if(
      _stereocenterMap.count(i) == 0  // not already in the map
      && molecule.getNumAdjacencies(i) > 1 // non-terminal
    ) {
      auto localRanking = molecule.rankPriority(i);

      _stereocenterMap[i] = std::make_shared<Stereocenters::CNStereocenter>(
        molecule.getGraph(),
        molecule.determineLocalGeometry(i, localRanking),
        i,
        localRanking
      );

      /* New stereocenters encountered at this point can have multiple
       * assignments, since some types of stereocenters are flatly ignored by
       * the candidate functions from Molecule, such as trigonal pyramidal
       * nitrogens. These are found here, though, and MUST be chosen randomly
       * according to the relative weights
       */
      std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
        _stereocenterMap[i]
      ) -> assignRandom();
    }
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
          cycleInternalAngles.at(angleCentralIndex - 1),
          angleAbsoluteVariance * looseningMultiplier
        );
      }

      // There's a missing triple with the previous approach, which is always
      setAngleBoundsIfEmpty(
        {{
          indexSequence.at(indexSequence.size() - 2),
          indexSequence.at(0),
          indexSequence.at(1)
        }},
        cycleInternalAngles.back(),
        angleAbsoluteVariance * looseningMultiplier
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
    if(smallestCycleMap.count(i) == 1) {
      if(smallestCycleMap.at(i) == 3) {
        return 2.5;
      } else if(smallestCycleMap.at(i) == 4) {
        return 1.7;
      } else if(smallestCycleMap.at(i) == 5) {
        return 1.3;
      }
    }

    return 1.0;
  };

  // Get all 1-3 and 1-4 information possible from stereocenters
  for(const auto& mapIterPair : _stereocenterMap) {
    const auto& stereocenterPtr = mapIterPair.second;

    stereocenterPtr->setModelInformation(
      *this,
      cycleMultiplierForIndex,
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
  for(const auto& stereocenterPtr : _molecule.getStereocenterList()) {
    if(stereocenterPtr -> type() == Stereocenters::Type::CNStereocenter) {
      auto cnPtr = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
        stereocenterPtr
      );

      AtomIndexType i = cnPtr->involvedAtoms().front();

      if(
        cnPtr -> getSymmetry() == Symmetry::Name::Tetrahedral
        && cycleData.numCycleFamilies(i) == 2
      ) {
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
                if(cycleOneVertices.count(iAdjacent)) {
                  firstAdjacents.push_back(iAdjacent);
                }

                if(cycleTwoVertices.count(iAdjacent)) {
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
    }
  }

  // Add default angles and dihedrals for all adjacent sequences
  addDefaultAngles();
  addDefaultDihedrals();
}

void MoleculeSpatialModel::setBondBoundsIfEmpty(
  const std::array<AtomIndexType, 2>& bondIndices,
  const double centralValue
) {
  double relativeVariance = bondRelativeVariance * _looseningMultiplier;

  /* Ensure correct use, variance for bounds should at MOST smaller than half
   * of the central value
   */
  assert(relativeVariance < 0.5 * centralValue);

  // C++17: try_emplace
  auto indexSequence = orderedIndexSequence<2>(bondIndices);
  if(_bondBounds.count(indexSequence) == 0) {
    _bondBounds.emplace(
      std::make_pair(
        indexSequence,
        ValueBounds {
          (1 - relativeVariance) * centralValue,
          (1 + relativeVariance) * centralValue
        }
      )
    );
  }
}

void MoleculeSpatialModel::setBondBoundsIfEmpty(
  const std::array<AtomIndexType, 2>& bondIndices,
  const ValueBounds& bounds
) {
  // C++17: try_emplace
  auto indexSequence = orderedIndexSequence<2>(bondIndices);
  if(_bondBounds.count(indexSequence) == 0) {
    _bondBounds.emplace(
      std::make_pair(
        indexSequence,
        bounds
      )
    );
  }
}

/*!
 * Adds the angle bounds to the model, but only if the information for that
 * set of indices does not exist yet.
 */
void MoleculeSpatialModel::setAngleBoundsIfEmpty(
  const std::array<AtomIndexType, 3>& angleIndices,
  const double centralValue,
  const double absoluteVariance
) {
  auto orderedIndices = orderedIndexSequence<3>(angleIndices);

  // C++17: try_emplace
  if(_angleBounds.count(orderedIndices) == 0) {
    _angleBounds.emplace(
      std::make_pair(
        orderedIndices,
        ValueBounds {
          StdlibTypeAlgorithms::clamp(
            centralValue - absoluteVariance,
            0.0,
            M_PI
          ),
          StdlibTypeAlgorithms::clamp(
            centralValue + absoluteVariance,
            0.0,
            M_PI
          )
        }
      )
    );
  }
}

void MoleculeSpatialModel::setAngleBoundsIfEmpty(
  const std::array<AtomIndexType, 3>& angleIndices,
  const ValueBounds& bounds
) {
  auto orderedIndices = orderedIndexSequence<3>(angleIndices);

  // C++17: try_emplace
  if(_angleBounds.count(orderedIndices) == 0) {
    _angleBounds.emplace(
      std::make_pair(
        orderedIndices,
        ValueBounds {
          StdlibTypeAlgorithms::clamp(
            bounds.lower,
            0.0,
            M_PI
          ),
          StdlibTypeAlgorithms::clamp(
            bounds.upper,
            0.0,
            M_PI
          )
        }
      )
    );
  }
}

/*!
 * Adds the dihedral bounds to the model, but only if the information for that
 * set of indices does not exist yet.
 */
void MoleculeSpatialModel::setDihedralBoundsIfEmpty(
  const std::array<AtomIndexType, 4>& dihedralIndices,
  const double lower,
  const double upper
) {
  assert(lower < upper);

  auto orderedIndices = orderedIndexSequence<4>(dihedralIndices);

  // C++17: try_emplace
  if(_dihedralBounds.count(orderedIndices) == 0) {
    _dihedralBounds.emplace(
      std::make_pair(
        orderedIndices,
        ValueBounds {
          StdlibTypeAlgorithms::clamp(
            lower,
            -M_PI,
            M_PI
          ),
          StdlibTypeAlgorithms::clamp(
            upper,
            -M_PI,
            M_PI
          )
        }
      )
    );
  }
}

void MoleculeSpatialModel::addDefaultAngles() {
  const AtomIndexType N = _molecule.numAtoms();
  for(AtomIndexType center = 0; center < N; ++center) {
    temple::forAllPairs(
      _molecule.getAdjacencies(center),
      [&](const AtomIndexType i, const AtomIndexType j) -> void {
        if(i != j) {
          setAngleBoundsIfEmpty(
            {{i, center, j}},
            0,
            M_PI
          );
        }
      }
    );
  }
}

void MoleculeSpatialModel::addDefaultDihedrals() {
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
            0,
            M_PI
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

DistanceBoundsMatrix MoleculeSpatialModel::makeBounds() const {
  DistanceBoundsMatrix bounds {_molecule.numAtoms()};

  for(const auto& bondPair : _bondBounds) {
    const auto& indices = bondPair.first;
    const auto& bondBounds = bondPair.second;

    assert(indices.front() != indices.back());

    /* Any elements in the bond bounds MUST improve the existing bounds and
     * may not cause contradictions.
     */
    assert(
      bounds.setLowerBound(indices.front(), indices.back(), bondBounds.lower)
      && bounds.setUpperBound(indices.front(), indices.back(), bondBounds.upper)
    );
  }

  for(const auto& anglePair : _angleBounds) {
    const auto& indices = anglePair.first;
    const auto& angleBounds = anglePair.second;

    assert(bondInformationIsPresent(bounds, indices));

    bounds.setLowerBound(
      indices.front(),
      indices.back(),
      CommonTrig::lawOfCosines(
        bounds.lowerBound(
          indices.front(),
          indices.at(1)
        ),
        bounds.lowerBound(
          indices.at(1),
          indices.back()
        ),
        angleBounds.lower
      )
    );

    bounds.setUpperBound(
      indices.front(),
      indices.back(),
      CommonTrig::lawOfCosines(
        bounds.upperBound(
          indices.front(),
          indices.at(1)
        ),
        bounds.upperBound(
          indices.at(1),
          indices.back()
        ),
        angleBounds.upper
      )
    );
  }

  for(const auto& dihedralPair : _dihedralBounds) {
    const auto& indices = dihedralPair.first;
    const auto& dihedralBounds = dihedralPair.second;

    assert(bondInformationIsPresent(bounds, indices));

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

    bounds.setLowerBound(
      indices.front(),
      indices.back(),
      CommonTrig::dihedralLength(
        bounds.lowerBound(
          indices.front(),
          indices.at(1)
        ),
        bounds.lowerBound(
          indices.at(1),
          indices.at(2)
        ),
        bounds.lowerBound(
          indices.at(2),
          indices.back()
        ),
        abAngleBounds.lower,
        bcAngleBounds.lower,
        dihedralBounds.lower // cis dihedral
      )
    );

    bounds.setUpperBound(
      indices.front(),
      indices.back(),
      CommonTrig::dihedralLength(
        bounds.upperBound(
          indices.front(),
          indices.at(1)
        ),
        bounds.upperBound(
          indices.at(1),
          indices.at(2)
        ),
        bounds.upperBound(
          indices.at(2),
          indices.back()
        ),
        abAngleBounds.upper,
        bcAngleBounds.upper,
        dihedralBounds.upper // trans dihedral
      )
    );
  }

  return bounds;
}

boost::optional<ValueBounds> MoleculeSpatialModel::coneAngle(
  const std::vector<AtomIndexType>& ligandIndices,
  const ValueBounds& coneHeightBounds
) const {
  return coneAngle(
    ligandIndices,
    coneHeightBounds,
    bondRelativeVariance * _looseningMultiplier,
    _molecule.getGraph()
  );
}

ValueBounds MoleculeSpatialModel::ligandDistance(
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

std::vector<DistanceGeometry::ChiralityConstraint> MoleculeSpatialModel::getChiralityConstraints() const {
  std::vector<DistanceGeometry::ChiralityConstraint> constraints;

  for(const auto& iterPair : _stereocenterMap) {
    const auto& stereocenterPtr = iterPair.second;

    auto chiralityPrototypes = stereocenterPtr -> chiralityConstraints();
    std::move(
      chiralityPrototypes.begin(),
      chiralityPrototypes.end(),
      std::back_inserter(constraints)
    );
  }

  return constraints;
}

void MoleculeSpatialModel::dumpDebugInfo() const {
  auto& logRef = Log::log(Log::Level::Debug);
  logRef << "MoleculeSpatialModel debug info" << std::endl;

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

struct MoleculeSpatialModel::ModelGraphWriter {
  /* Settings to determine appearance */
  // Color maps
  const std::map<
    std::string,
    std::string
  > elementBGColorMap {
    {"H", "white"},
    {"C", "gray"},
    {"N", "blue"},
    {"O", "red"}
  };

  const std::map<
    std::string,
    std::string
  > elementTextColorMap {
    {"H", "black"},
    {"C", "white"},
    {"N", "white"},
    {"O", "white"}
  };

  const std::map<
    BondType,
    std::string
  > bondTypeDisplayString {
    {BondType::Single, R"(color = "black")"},
    {BondType::Double, R"(color = "black:invis:black")"},
    {BondType::Triple, R"(color = "black:invis:black:invis:black")"},
    {BondType::Quadruple, R"(label = "4")"},
    {BondType::Quintuple, R"(label = "5")"},
    {BondType::Sextuple, R"(label = "6")"},
    {BondType::Aromatic, R"(style = "dashed")"},
    {BondType::Eta, R"(style = "dotted")"}
  };

  /* State */
  // We promise to be good and not change anything
  const GraphType* const graphPtr;
  const MoleculeSpatialModel& spatialModel;

/* Constructor */
  ModelGraphWriter(
    const GraphType* passGraphPtr,
    const MoleculeSpatialModel& spatialModel
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
    os << "graph [fontname = \"Arial\", layout = neato];\n"
      << "node [fontname = \"Arial\", shape = circle, style = filled];\n"
      << "edge [fontname = \"Arial\"];\n";

    /* Additional nodes */
    for(const auto& stereocenterIterPair : spatialModel._stereocenterMap) {
      const auto& mappingIndex = stereocenterIterPair.first;
      const auto& stereocenterPtr = stereocenterIterPair.second;

      if(
        stereocenterPtr->type() == Stereocenters::Type::EZStereocenter
      ) {
        // Avoid writing EZStereocenters twice
        if(*stereocenterPtr->involvedAtoms().rbegin() == mappingIndex) {
          continue;
        }

        char ezState = 'u';
        if(stereocenterPtr -> assigned()) {
          if(stereocenterPtr -> assigned().value() == 1) {
            ezState = 'Z';
          } else {
            ezState = 'E';
          }
        }

        os << "EZ" << temple::condenseIterable(
            stereocenterPtr -> involvedAtoms(),
            ""
          ) << R"( [label=")" << ezState
          << R"(", fillcolor="tomato", shape="square", fontcolor="white", )"
          << R"(tooltip=")"
          << stereocenterPtr -> info()
          << R"("];)" << "\n";
      } else {
        const auto cnPtr = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
          stereocenterPtr
        );

        os << "CN" << temple::condenseIterable(
          stereocenterPtr -> involvedAtoms(),
          ""
        ) << R"( [label=")" << Symmetry::name(cnPtr -> getSymmetry())
          << R"(", fillcolor="steelblue", shape="square", fontcolor="white", )"
          << R"(tooltip=")" << cnPtr -> info()
          << R"("];)" << "\n";
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
    if(elementBGColorMap.count(symbolString) != 0u) {
      os << R"(, fillcolor=")" << elementBGColorMap.at(symbolString) << R"(")";
    } else { // default
      os << R"(, fillcolor="white")";
    }
    if(elementTextColorMap.count(symbolString) != 0u) {
      os << R"(, fontcolor=")" << elementTextColorMap.at(symbolString) << R"(")";
    } else { // default
      os << R"(, fontcolor="orange")";
    }

    // Font sizing
    if(symbolString == "H") {
      os << ", fontsize=10, width=.3, fixedsize=true";
    }

    // Any angles this atom is the central atom in
    std::vector<std::string> angleStrings;
    for(const auto& angleIterPair : spatialModel._angleBounds) {
      const auto& indexSequence = angleIterPair.first;
      const auto& angleBounds = angleIterPair.second;

      if(indexSequence.at(1) == vertexIndex) {
        angleStrings.emplace_back(
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

    if(angleStrings.empty()) {
      os << R"(, tooltip="no angles here")";
    } else {
      os << R"(, tooltip="angles: )" << temple::condenseIterable(
        angleStrings,
        "&#10;"s
      ) << R"(")";
    }

    os << "]\n";

    // Are there any additional vertices we ought to connect to?
    if(spatialModel._stereocenterMap.count(vertexIndex) == 1) {
      const auto& stereocenterPtr = spatialModel._stereocenterMap.at(vertexIndex);
      if(stereocenterPtr -> type() == Stereocenters::Type::CNStereocenter) {
        os << ";CN";
      } else {
        os << ";EZ";
      }

      os << temple::condenseIterable(stereocenterPtr -> involvedAtoms(), "")
        << " -- " << vertexIndex << R"([color="gray", dir="forward", len="2"])";
    }
  }

  // Edge options
  void operator() (std::ostream& os, const EdgeIndexType& edgeIndex) const {
    const AtomIndexType source = boost::source(edgeIndex, *graphPtr);
    const AtomIndexType target = boost::target(edgeIndex, *graphPtr);

    os << "[";

    // Bond Type display options
    auto bondType = (*graphPtr)[edgeIndex].bondType;
    if(bondTypeDisplayString.count(bondType) != 0u) {
      os << bondTypeDisplayString.at(bondType);
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

void MoleculeSpatialModel::writeGraphviz(const std::string& filename) const {
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
boost::optional<ValueBounds> MoleculeSpatialModel::coneAngle(
  const std::vector<AtomIndexType>& baseConstituents,
  const ValueBounds& coneHeightBounds,
  const double bondRelativeVariance,
  const GraphType& graph
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
  /* TODO This works, but it might be preferable to just iterate through
   * combinations of RCs instead, the full cycle MUST be there. Unsure if this
   * may be prohibitively expensive, though.
   */
  // For this cycle interpretation, we need to ignore eta bonds
  Cycles cycles {graph, true};

  /*for(const auto cyclePtr : cycles) {
    std::cout << "Cycle " << temple::condenseIterable(
      makeRingIndexSequence(
        Cycles::edgeVertices(cyclePtr)
      )
    ) << std::endl;
  }*/

  /* This is essentially an if - only one cycle can consist of exactly the base
   * constituents
   */
  for(
    const auto cyclePtr :
    cycles.iterate(Cycles::predicates::ConsistsOf {baseConstituents})
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
   * For some types of ligands, we can still figure out a cone angle. If the
   * ligand group is actually a path in which any intermediate atom merely
   * connects the subsequent atoms (i.e. the path is not branched, there are no
   * cycles), and if all the involved intermediate geometries constituting the
   * longest path have only one distinct angle value, then we can create a
   * conformational model anyway.
   */

  // Try to find a path including all atoms without any cross-bonds
  // Wrap baseConstituents with an TinyUnorderedSet for easier member searches
  temple::TinyUnorderedSet<AtomIndexType> ligandIndicesSet;
  ligandIndicesSet.data = baseConstituents;

  std::deque<AtomIndexType> pathSequence;
  pathSequence.push_back(baseConstituents.front());

  auto adjacentIters = boost::adjacent_vertices(baseConstituents.front(), graph);
  std::vector<AtomIndexType> initialAdjacents;
  std::copy(
    adjacentIters.first,
    adjacentIters.second,
    std::back_inserter(initialAdjacents)
  );

  temple::TinyUnorderedSet<AtomIndexType> initialAdjacentsSet;
  initialAdjacentsSet.data = initialAdjacents;

  auto intersection = temple::unorderedSetIntersection(
    ligandIndicesSet,
    initialAdjacentsSet
  );

  // Early exit, there can be no path involving all without cross-bonds
  if(intersection.size() > 2 || intersection.size() == 0) {
    return boost::none;
  }

  bool frontFinished = false;
  bool backFinished = false;
  bool crossEdgeDiscovered = false;

  if(intersection.size() == 1) {
    pathSequence.push_front(intersection.data.front());
    backFinished = true;
  }

  if(intersection.size() == 2) {
    pathSequence.push_front(intersection.data.front());
    pathSequence.push_back(intersection.data.back());
  }

  while(pathSequence.size() != ligandIndicesSet.size()) {
    if(frontFinished && backFinished) {
      break;
    }

    // Try to expand at front
    if(!frontFinished) {
      frontFinished = true;
      unsigned count = 0;
      for(
        const AtomIndexType adjacent :
        RangeForTemporary<GraphType::adjacency_iterator>(
          boost::adjacent_vertices(pathSequence.front(), graph)
        )
      ) {
        if(ligandIndicesSet.count(adjacent)) {
          if(adjacent != pathSequence.at(1)) {
            pathSequence.push_front(adjacent);
            frontFinished = false;
          }

          ++count;
        }
      }

      if(count > 2) {
        crossEdgeDiscovered = true;
        break;
      }
    }

    if(!backFinished) {
      backFinished = true;
      unsigned count = 0;
      for(
        const auto adjacent :
        RangeForTemporary<GraphType::adjacency_iterator>(
          boost::adjacent_vertices(pathSequence.back(), graph)
        )
      ) {
        if(ligandIndicesSet.count(adjacent)) {
          if(adjacent != pathSequence.at(pathSequence.size() - 2)) {
            pathSequence.push_back(adjacent);
            backFinished = false;
          }

          ++count;
        }
      }

      if(count > 2) {
        crossEdgeDiscovered = true;
        break;
      }
    }
  }

  if(pathSequence.size() == ligandIndicesSet.size() && !crossEdgeDiscovered) {
    // TODO model the path and determine a cone angle!
    throw std::logic_error("Modeling of path ligands is not implemented.");
  }

  return boost::none;
}

double MoleculeSpatialModel::spiroCrossAngle(const double alpha, const double beta) {
  // The source of this equation is explained in documents/
  return std::acos(
    -1 * std::cos(alpha / 2) * std::cos(beta / 2)
  );
}

ValueBounds MoleculeSpatialModel::ligandDistanceFromCenter(
  const std::vector<AtomIndexType>& ligandIndices,
  const AtomIndexType centralIndex,
  const double bondRelativeVariance,
  const GraphType& graph
) {
  assert(ligandIndices.size() > 0);

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

} // namespace DistanceGeometry

} // namespace molassembler
