#include "Molecule.h"

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
 */

namespace MoleculeManip {

/* Constructors */
Molecule::Molecule(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) {
  // store ElementTypes
  _elements.push_back(a);
  _elements.push_back(b);
  // update _adjacencies
  _adjacencies.addSlot();
  _adjacencies.addSlot();
  _adjacencies.addAdjacency(0, 1);
  // update _edges
  _edges.add(Edge(
    0,
    1,
    bondType
  ));
}

Molecule::Molecule(
  const Delib::ElementTypeCollection& elements,
  const AdjacencyList& adjacencies,
  const EdgeList& edges
) : 
  _elements(elements),
  _adjacencies(adjacencies),
  _edges(edges)
{}

Molecule::Molecule(
  const Delib::ElementTypeCollection& elements,
  const Delib::PositionCollection& positions,
  const AdjacencyList& adjacencies,
  const EdgeList& edges
) : 
  _elements(elements),
  _positions(positions),
  _adjacencies(adjacencies),
  _edges(edges)
{}

/* Private member functions */
bool Molecule::_validAtomIndex(const AtomIndexType& a) const {
  return (
    a < _adjacencies.size()
    && _elements[a] != Delib::ElementType::none
  );
}

bool Molecule::_validAtomIndices(
  const AtomIndexType& a,
  const AtomIndexType& b
) const {
  return (
    _validAtomIndex(a)
    && _validAtomIndex(b)
    && a < b
  );
}

AtomIndexType Molecule::addAtom(
  const Delib::ElementType& elementType,
  const AtomIndexType& bondedToIndex,
  const BondType& bondType
) {
  auto addedIndex = _adjacencies.addSlot();

  assert(bondedToIndex < addedIndex);
  _adjacencies.addAdjacency(
    bondedToIndex,
    addedIndex
  );

  _edges.add(Edge(
    bondedToIndex,
    addedIndex,
    bondType
  ));

  _elements.push_back(elementType);

  return addedIndex;
}

void Molecule::addBond(
  const AtomIndexType& a,
  const AtomIndexType& b,
  const BondType& bondType
) {
  assert(_validAtomIndices(a, b));
  _adjacencies.addAdjacency(
    a,
    b
  );

  _edges.add(Edge(
    a,
    b,
    bondType
  ));
}

void Molecule::removeAtom(const AtomIndexType& a) {
  assert(_validAtomIndex(a));

  // must update edges and adjacencies 
  auto bonded_to = _adjacencies.getAdjacencies(a);

  // erase all edges to and from this atom
  for(const auto& bondedAtomIndex : bonded_to) {
    _edges.remove(
      std::min(a, bondedAtomIndex),
      std::max(a, bondedAtomIndex)
    );
  }

  // remove all other mentions in _adjacencies
  for(const auto& bondedAtomIndex : bonded_to) {
    _adjacencies.removeAdjacency(a, bondedAtomIndex);
  }

  // set _atom_exists to false for this index
  _elements[a] = Delib::ElementType::none;
}

void Molecule::removeBond(
  const AtomIndexType& a,
  const AtomIndexType& b
) {
  assert(_validAtomIndices(a, b));
  _adjacencies.removeAdjacency(a, b);

  _edges.remove(
    a,
    b
  );
}

Delib::ElementType Molecule::getElementType(
  const AtomIndexType& a
) const {
  assert(_validAtomIndex(a));
  return _elements[a];
}

std::vector<AtomIndexType> Molecule::getBondedAtomIndices(
  const AtomIndexType& a
) const {
  assert(_validAtomIndex(a));
  return _adjacencies.getAdjacencies(a);
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
> Molecule::rankCIPPriority(
  const AtomIndexType& a,
  const std::vector<AtomIndexType>& excludeAdjacent 
) const {
  auto toRank = _adjacencies.getAdjacencies(a);
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
    std::set<AtomIndexType> visitedSet,
    std::vector<AtomIndexType> seeds,
    std::multiset<int, std::greater<int> > Zs
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
  
  // sort toRank according to CIP
  std::sort(
    toRank.begin(),
    toRank.end(),
    [&a, this, getZ, BFSIterate, &equalPairs](
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
          lhsSeeds.size() != 0
          || rhsSeeds.size() != 0
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

std::ostream& operator << (
  std::ostream& os,
  const Molecule& mol
) {
  os << mol._adjacencies.size() << " adjacencies, "
    << mol._edges.size() << " edges." << std::endl;
  return os;
}

} // eo namespace
