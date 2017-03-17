#include "Molecule.h"
#include "CNStereocenter.h"
#include "EZStereocenter.h"
#include "BondDistance.h"
#include "symmetry_information/Symmetries.h"
#include "CommonTrig.h"

#include "GraphDistanceMatrix.h"
#include "AdjacencyListAlgorithms.h"
#include "AdjacencyMatrix.h"

#include "DistanceGeometry/DistanceBoundsMatrix.h"

// Delib
#include "ElementInfo.h"

/* TODO
 * - implement remaining functions, importantly the operator ==
 * - consider how to interface GraphFeatures with Molecule, there is a need for
 *   information flow between both:: Molecule operator == needs to be able
 *   to compare GraphFeatures and GraphFeatures need access to element types
 *   and edges
 * - Molecule operator == can be implemented in myriad ways, suggest some 
 *   heuristics before going for testing graph isomorphism:
 *
 *   1. #vertices, #edges O(1)
 *   2. comparison of GraphFeatures 
 *      - same amount
 *      - exactly the same
 *   3. isomorphism
 *
 *   maybe also consider the Wiener index, global clustering and algebraic
 *   connectivity
 *
 *   (in progress)
 * - alter all functions to respect the internal-external divide with internal 
 *   indices always being contiguous!
 *
 */

namespace MoleculeManip {

using namespace DistanceGeometry;

/* Constructors --------------------------------------------------------------*/
Molecule::Molecule(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) {
  // update _adjacencies
  _adjacencies.addAtom(a);
  _adjacencies.addAtom(b);
  _adjacencies.addBond(0, 1, bondType);
}

Molecule::Molecule(const AdjacencyList& adjacencies) 
: _adjacencies(adjacencies)
{ 
  _detectStereocenters();
}

/* Private member functions --------------------------------------------------*/
// TODO deprecate and remove as soon as BFSConstraintCollector is ready
std::vector<DistanceConstraint> Molecule::_createConstraint(
  const std::vector<AtomIndexType>& chain
) const {
  auto& i = chain.front();
  auto& j = chain.back();

  switch(chain.size()) {
    case 2: {
      auto distance = Bond::calculateBondDistance(
        getElementType(i),
        getElementType(j),
        getBondType(i, j)
      );
      return {
        DistanceConstraint(
          i,
          j,
          0.99 * distance,
          1.01 * distance
        )
      };
    }
    case 3: {
      auto& intermediate = chain[1];
      // TODO continue here
      // no considerations of charge at all so far, not even VSEPR
      boost::optional<double> angle;
      // determine which symmetry applies at intermediate position, angle
      auto nLigands = getBondedAtomIndices(intermediate).size();
      if(nLigands == 2) {
        if(getElementType(intermediate) == Delib::ElementType::O) {
          angle = 104.5;
        } else {
          angle = Symmetry::angleFunction(Symmetry::Name::Linear)(0, 1);
        }
      } else if(nLigands == 3) {
        angle = Symmetry::angleFunction(Symmetry::Name::TrigonalPlanar)(0, 1);
      } else if(nLigands == 4) {
        angle = Symmetry::angleFunction(Symmetry::Name::Tetrahedral)(0, 1);
      } 

      if((bool) angle) {
        auto a = Bond::calculateBondDistance(
          getElementType(i),
          getElementType(intermediate),
          getBondType(i, intermediate)
        );
        auto b = Bond::calculateBondDistance(
          getElementType(intermediate),
          getElementType(j),
          getBondType(intermediate, j)
        );
        auto distance = CommonTrig::lawOfCosines(
          a,
          b,
          angle.value() * M_PI / 180.0
        );

        return {
          DistanceConstraint(
            i,
            j,
            0.95 * distance,
            1.05 * distance
          )
        };
      } else {
        return {};
      }
    }
    case 4: {
      // TODO implement
      return {};
    }
    default: {
      return {
        DistanceConstraint(
          i,
          j,
          (
            AtomInfo::vdwRadius(getElementType(i))
            + AtomInfo::vdwRadius(getElementType(j))
          ),
          100
        )
      };
    }
  }
}

void Molecule::_detectStereocenters() {
  // Find CNStereocenters
  for(unsigned i = 0; i < _adjacencies.size(); i++) {
    if(
      /* TODO this is no longer a valid way of checking how many ligands there are
       * -> eta bonds exist!
       */
      _adjacencies.getNumAdjacencies(i) >= 3 
    ) {
      // TODO geometry determination for transition metal complexes?

      // Determine the local geometry
      auto localGeometryName = _determineLocalGeometry(i);
      auto rankResultPair = rankPriority(i);

      // Construct a Stereocenter here
      std::shared_ptr<
        Stereocenters::CNStereocenter
      > newStereocenter = std::make_shared<
        Stereocenters::CNStereocenter
      >(
        localGeometryName,
        i,
        rankResultPair.first,
        rankResultPair.second
      );

      if(newStereocenter -> assignments() > 1) {
        _stereocenters.add(
          std::move(newStereocenter)
        );
      }
    }
  }

  /* TODO
   * - Will need refinement to not instantiate EZStereocenters in small cycles
   *   (up to a preset size, maybe around 8 or so?)
   */
  // Find EZStereocenters
  const auto& graphRef = _adjacencies.access();
  for(
    const auto& edgeIndex : 
    RangeForTemporary<GraphType::edge_iterator>(
      boost::edges(graphRef)
    )
  ) {
    auto source = boost::source(edgeIndex, graphRef),
         target = boost::target(edgeIndex, graphRef);

    if(
      graphRef[edgeIndex].bondType == BondType::Double
      && _adjacencies.getNumNonEtaAdjacencies(source) == 3
      && _adjacencies.getNumNonEtaAdjacencies(target) == 3
    ) {
      // Calculate Priorities for each's substituents
      auto sourceSubstituentsRanking = rankPriority(
        source,
        {target} // exclude edge sharing neighbor
      );

      // If the source's substituents are unequal (no equal pair sets)
      if(sourceSubstituentsRanking.second.size() == 0) {
        auto targetSubstituentsRanking = rankPriority(
          target,
          {source} // exclude edge sharing neighbor
        );

        // target must also have no equal pairs
        if(targetSubstituentsRanking.second.size() == 0) {
          // Instantiate an EZStereocenter there!
          _stereocenters.add(
            std::make_shared<
              Stereocenters::EZStereocenter
            >(
              source,
              sourceSubstituentsRanking.first,
              target,
              targetSubstituentsRanking.first
            )
          );
        }
      }
    }
  }

}

std::vector<LocalGeometry::LigandType> Molecule::_reduceToLigandTypes(
  const AtomIndexType& index
) const {
  /* TODO 
   * - No L, X determination. Although, will L, X even be needed for metals?
   *   Maybe only for OZ and NVE determination...
   */
  /* VSEPR formulation is that geometry is a function of 
   * - localized charge of central atom
   * - atom type of central atom, neighbors
   * - bond types to neighbors
   */

  // call this only on non-terminal atoms
  assert(_adjacencies.getNumAdjacencies(index) > 1);

  // first basic stuff for VSEPR, later L and X for transition metals
  // geometry inference does not care if the substituents are somehow 
  // connected (unless in later models the entire structure is considered)
  std::vector<LocalGeometry::LigandType> ligands;

  for(const auto& adjacentIndex: _adjacencies.iterateAdjacencies(index)) {
    ligands.push_back(
      LocalGeometry::LigandType {
        0, 0, {{  // L and X are 0 since only VSEPR is considered for now
          _adjacencies.getElementType(adjacentIndex),
          getBondType(index, adjacentIndex)
        }}
      }
    );
  }

  return ligands;
}

Symmetry::Name Molecule::_determineLocalGeometry(
  const AtomIndexType& index
) const {
  assert( _adjacencies.getNumAdjacencies(index) > 1); 

  auto ligandsVector = _reduceToLigandTypes(index);
  // TODO this below is invalid for metals!
  unsigned nSites = _adjacencies.getNumAdjacencies(index);
  int formalCharge = 0;

  return LocalGeometry::VSEPR::determineGeometry(
    _adjacencies.getElementType(index),
    nSites,
    ligandsVector,
    formalCharge
  );
}

std::map<AtomIndexType, Symmetry::Name> Molecule::_determineLocalGeometries() const {
  std::map<AtomIndexType, Symmetry::Name> symmetryMap;

  for(AtomIndexType i = 0; i < _adjacencies.size(); i++) {
    if(_adjacencies.getNumAdjacencies(i) > 1) {
      symmetryMap[i] = _determineLocalGeometry(i);
    }
  }

  return symmetryMap;
}

/* Public modifiers ----------------------------------------------------------*/
AtomIndexType Molecule::addAtom(
  const Delib::ElementType& elementType,
  const AtomIndexType& bondedToIndex,
  const BondType& bondType
) {
  auto addedIndex = _adjacencies.addAtom(elementType);

  _adjacencies.addBond(
    bondedToIndex,
    addedIndex,
    bondType
  );

  return addedIndex;
}

void Molecule::addBond(
  const AtomIndexType& a,
  const AtomIndexType& b,
  const BondType& bondType
) {
  _adjacencies.addBond(a, b, bondType);
}

void Molecule::removeAtom(const AtomIndexType& a) {
  assert(_validAtomIndex(a));

  // copy out the adjacency list
  auto adjacencyListCopy = _adjacencies;

  // Remove all vertices to and from the atom, plus the atom itself
  adjacencyListCopy.removeAtom(a);

  // is this still a connected molecule or have we split it in two?
  if(AdjacencyListAlgorithms::numConnectedComponents(adjacencyListCopy) != 1) {
    throw std::logic_error(
      "The selected atom removal would lead to a molecule split!"
    );
  }

  // if nothrow, overwrite internal with modified adjacencyList
  _adjacencies = adjacencyListCopy;

  // TODO URGENT StereocenterList index invalidation update
  // just kinda useless if Stereocenter is going to change a lot
  // additionally, if an atom is removed, maybe ranking changes? some update 
  // algorithm is needed
}

void Molecule::removeBond(
  const AtomIndexType& a,
  const AtomIndexType& b
) {
  assert(_validAtomIndex(a) && _validAtomIndex(b));

  auto adjacencyListCopy = _adjacencies;

  adjacencyListCopy.removeBond(a, b);

  if(AdjacencyListAlgorithms::numConnectedComponents(adjacencyListCopy) != 1) {
    throw std::logic_error(
      "The selected bond removal would lead to a molecule split!"
    );
  }

  _adjacencies = adjacencyListCopy;

  // No minimization necessary here -> The removed bond does not lead to a 
  // disconnected atom (otherwise we would have thrown a few lines back). But a
  // bond removal can require a stereocenter update => TODO
}

/* Public information --------------------------------------------------------*/
void Molecule::dumpGraphviz(const std::string& filename) const {
  _adjacencies.dumpGraphviz(filename);
}

bool Molecule::_validAtomIndex(const AtomIndexType& a) const {
  return (
    a < _adjacencies.size()
  );
}

const AdjacencyList& Molecule::getAdjacencyList() const {
  return _adjacencies;
}

std::vector<AtomIndexType> Molecule::getBondedAtomIndices(
  const AtomIndexType& a
) const {
  assert(_validAtomIndex(a));
  return _adjacencies.getAdjacencies(a);
}

BondType Molecule::getBondType(
  const AtomIndexType& a,
  const AtomIndexType& b
) const {
  auto edgeOption = _adjacencies.getBondType(a, b);

  if(!edgeOption) {
    throw std::logic_error("The atoms requested are not bonded!");
  }

  return edgeOption.value();
}

DistanceGeometry::DistanceBoundsMatrix Molecule::getDistanceBoundsMatrix() const {

  const double oneTwoDelta = 0.01;

  DistanceGeometry::DistanceBoundsMatrix distanceBounds(getNumAtoms());

  /* enter all 1-2 distance bounds into the matrix for reference, they will be
   * required often
   */

  for(AtomIndexType i = 0; i < _adjacencies.size(); i++) {
    for(const auto& j: getBondedAtomIndices(i)) {
      // avoid duplicate work by avoiding i > j cases
      if(i > j) continue;

      // enter the constraints for i, j
      const auto distance = Bond::calculateBondDistance(
        getElementType(i),
        getElementType(j),
        getBondType(i, j)
      );

      distanceBounds.upperBound(i, j) = (1 + oneTwoDelta) * distance;
      distanceBounds.lowerBound(i, j) = (1 - oneTwoDelta) * distance;
    }
  }

  return distanceBounds;
}

std::vector<AdjacencyList::ExplicitEdge> Molecule::getEdges() const {
  return _adjacencies.getEdges();
}

Delib::ElementType Molecule::getElementType(
  const AtomIndexType& a
) const {
  return _adjacencies.getElementType(a);
}

unsigned Molecule::getNumAtoms() const {
  return _adjacencies.nAtoms();
}

unsigned Molecule::getNumBonds() const {
  return _adjacencies.nBonds();
}

unsigned Molecule::hydrogenCount(const AtomIndexType& a) const {
  assert(_validAtomIndex(a));
  unsigned count = 0;

  for(const auto& adjacentIndex : _adjacencies[a]) {
    if(getElementType(adjacentIndex) == Delib::ElementType::H) {
      count += 1;
    }
  }

  return count;
}

/* TODO
 * - does not treat correctly:
 *   - cycles
 *   - stereocenters (Z over E, R over S (?))
 *   - double and triple bond ghost atom splitting
 * - unsure about sub-lists. is this approach even remotely correct?
 * - FUCK CIP rules -> maybe just use the unsigned values of assignments in
 *   GraphFeatures and rank branches with that.
 * - test
 */
std::pair<
  std::vector<AtomIndexType>, // the sorted list of substituent priorities
  std::set< // a set of pairs of AtomIndexTypes that are EQUAL
    std::pair<
      AtomIndexType,
      AtomIndexType
    >
  >
> Molecule::rankPriority(
  const AtomIndexType& a,
  const std::vector<AtomIndexType>& excludeAdjacent 
) const {
  auto toRank = getBondedAtomIndices(a);
  std::set<
    std::pair<
      AtomIndexType,
      AtomIndexType
    >
  > equalPairs;

  // remove excludes from toRank
  toRank.erase(
    std::remove_if(
      toRank.begin(),
      toRank.end(),
      [&excludeAdjacent](const AtomIndexType& atomIndex) {
        return std::find(
          excludeAdjacent.begin(),
          excludeAdjacent.end(),
          atomIndex
        ) != excludeAdjacent.end();
      }
    ),
    toRank.end()
  );

  auto getZ = [this](const AtomIndexType& a) -> int {
    return static_cast<int>(getElementType(a));
  };

  auto BFSIterate = [this, getZ](
    std::set<AtomIndexType>& visitedSet,
    std::vector<AtomIndexType>& seeds,
    std::multiset<int, std::greater<int> >& Zs
  ) {
    std::vector<AtomIndexType> newSeeds;
    for(const auto& index : seeds) {
      auto adjacent = getBondedAtomIndices(index);
      for(const auto& potentialSeed : adjacent) {
        // if not in visitedSet
        if(visitedSet.count(potentialSeed) == 0) {
          // add it to the set and new seeds
          visitedSet.insert(potentialSeed);
          newSeeds.push_back(potentialSeed);

          // add it's Z to Zs
          Zs.insert(getZ(potentialSeed));
        } // else skip
      }
    }

    // overwrite seeds
    seeds = newSeeds;
  };
  
  // sort toRank according to CIP-like rules
  std::sort(
    toRank.begin(),
    toRank.end(),
    [&a, this, &getZ, &BFSIterate, &equalPairs](
      const AtomIndexType& lhs,
      const AtomIndexType& rhs
    ) {
      std::set<AtomIndexType> lhsVisited = {a}, rhsVisited = {a};
      std::vector<AtomIndexType> lhsSeeds = {lhs}, rhsSeeds = {rhs};
      std::multiset<
        int,
        std::greater<int> // in CIP, list of Z is ordered DESC
      > lhsZ = { getZ(lhs) }, rhsZ = { getZ(rhs) };
      while(
          lhsSeeds.size() > 0
          || rhsSeeds.size() > 0
      ) {
        // compare lists
        if(lhsZ < rhsZ) return true;
        else if(lhsZ > rhsZ) return false;

        // iterate along the bonds
        BFSIterate(lhsVisited, lhsSeeds, lhsZ);
        BFSIterate(rhsVisited, rhsSeeds, rhsZ);
      }

      // all equal -> add to equalPairs
      equalPairs.emplace(
        std::min(lhs, rhs),
        std::max(lhs, rhs)
      );
      return false;
    }
  );

  return {
    toRank,
    equalPairs
  };
}

/* Operators -----------------------------------------------------------------*/
std::ostream& operator << (
  std::ostream& os,
  const Molecule& mol
) {
  if(mol._stereocenters.size() > 0) {
    os << "Stereocenter information:\n";
    for(const auto& stereocenterPtr: mol._stereocenters) {
      os << stereocenterPtr << std::endl;
    }
  }
  
  return os;
}

} // eo namespace
