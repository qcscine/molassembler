#include "Molecule.h"
#include "BondDistance.h"
#include "symmetry_information/Symmetries.h"
#include "CommonTrig.h"

#include "GraphDistanceMatrix.h"
#include "AdjacencyListAlgorithms.h"
#include "AdjacencyMatrix.h"

#include "TreeAlgorithms.h"

// Delib
#include "Delib/ElementInfo.h"

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
: _adjacencies(adjacencies), 
  stereocenters(adjacencies.detectStereocenters())
{}

Molecule::Molecule(
  const AdjacencyList& adjacencies,
  const StereocenterList& stereocenters
) : _adjacencies(adjacencies),
    stereocenters(stereocenters)
{}

/* Private member functions --------------------------------------------------*/

bool Molecule::_validAtomIndex(const AtomIndexType& a) const {
  return (
    a < _adjacencies.numAtoms()
  );
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

void Molecule::changeElementType(
  const AtomIndexType& a,
  const Delib::ElementType& elementType
) {
  _adjacencies.changeElementType(a, elementType);
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

void Molecule::updateStereocenters() {
  stereocenters = _adjacencies.detectStereocenters();
}

/* Public information --------------------------------------------------------*/
void Molecule::dumpGraphviz(const std::string& filename) const {
  _adjacencies.dumpGraphviz(filename);
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

std::vector<AdjacencyList::ExplicitEdge> Molecule::getEdges() const {
  return _adjacencies.getEdges();
}

Delib::ElementType Molecule::getElementType(
  const AtomIndexType& a
) const {
  return _adjacencies.getElementType(a);
}

unsigned Molecule::getNumAtoms() const {
  return _adjacencies.numAtoms();
}

unsigned Molecule::getNumBonds() const {
  return _adjacencies.numBonds();
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


/* Operators -----------------------------------------------------------------*/
std::ostream& operator << (
  std::ostream& os,
  const Molecule& mol
) {
  if(!mol.stereocenters.empty()) {
    os << "Stereocenter information:\n";
    for(const auto& stereocenterPtr: mol.stereocenters) {
      os << stereocenterPtr << std::endl;
    }
  }
  
  return os;
}

} // namespace MoleculeManip
