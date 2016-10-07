#include "Molecule.h"

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

} // eo namespace
