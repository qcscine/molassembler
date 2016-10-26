#ifndef INCLUDE_MOLECULE_H
#define INCLUDE_MOLECULE_H

// STL
#include <iostream>
#include <set>

// Delib
#include "Types/PositionCollection.h"
#include "Types/ElementTypeCollection.h"

// Custom headers
#include "AdjacencyList.h"
#include "EdgeList.h"
#include "StereocenterList.h"

namespace MoleculeManip {

class Molecule {
private:
/* private members */

  /* A Molecule conceptually contains a graph:
   * - Atoms are vertices (and thus have values)
   * - Bonds are edges (and thus weighted)
   * - The ensuing graph is
   *   - connected: a path from any node to any other exists
   *   - sparse: few edges present compared to the number of possible 
   *     edges -> use adjacency list instead of adjacency matrix
   */

  // The set of QC data on the atoms
  Delib::ElementTypeCollection _elements;
  Delib::PositionCollection _positions;
  // The information on interconnectedness of the atoms
  AdjacencyList _adjacencies;
  EdgeList _edges;
  StereocenterList _features;
  
  /* Private member functions */
  bool _validAtomIndex(const AtomIndexType& a) const;
  bool _validAtomIndices(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;

public:
/* Constructors */
  Molecule(
    const Delib::ElementType& a,
    const Delib::ElementType& b,
    const BondType& bondType
  ); 

  Molecule(
    const Delib::ElementTypeCollection& elements,
    const AdjacencyList& adjacencies,
    const EdgeList& edges
  );

  Molecule(
    const Delib::ElementTypeCollection& elements,
    const Delib::PositionCollection& positions,
    const AdjacencyList& adjacencies,
    const EdgeList& edges
  );

  AtomIndexType addAtom(
    const Delib::ElementType& elementType,
    const AtomIndexType& bondedToIndex,
    const BondType& bondType
  );

  void addBond(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const BondType& bondType
  ); 

  void removeAtom(const AtomIndexType& a); 

  void removeBond(
    const AtomIndexType& a,
    const AtomIndexType& b
  );

  /* Information retrieval */
  Delib::ElementType getElementType(
    const AtomIndexType& a
  ) const;
  /*bool bond_exists(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;
  bool bond_exists(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const BondType& bond_type
  ) const;
  BondType get_bond_type(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;
  std::vector<
    std::pair<
      AtomIndexType,
      BondType
    >
  > get_bond_pairs(const AtomIndexType& a) const;
  */

  std::vector<AtomIndexType> getBondedAtomIndices(
    const AtomIndexType& a
  ) const;

  std::pair<
    std::vector<AtomIndexType>, // the sorted list of substituent priorities
    std::set< // a set of pairs of AtomIndexTypes that are EQUAL
      std::pair<
        AtomIndexType,
        AtomIndexType
      >
    >
  > rankCIPPriority(
    const AtomIndexType& a,
    const std::vector<AtomIndexType>& excludeAdjacent
  ) const;

  /* Testing */
  std::pair<bool, std::string> validate() const noexcept;

/* Operators */
  /* An efficient implementation of the following two is imperative.
   * Some ideas for fast differentiation can probably be found from the 
   * wikipedia category "Graph invariants"
   */
  bool operator == (const Molecule& b) const;
  bool operator != (const Molecule& b) const;

/* Friends */
  /* Output stream operator for easier debugging. */
  friend std::ostream& operator<<(std::ostream& os, const Molecule& mol);
};

}

#endif
