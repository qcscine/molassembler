#ifndef INCLUDE_MOLECULE_H
#define INCLUDE_MOLECULE_H

// STL
#include <iostream>
#include <set>

// Custom headers
#include "AdjacencyList.h"
#include "StereocenterList.h"
#include "BondDistance.h"

namespace MoleculeManip {

/* TODO
 * - Figure out how to handle charge, where we need to pay attention and what 
 *   the minimal way is
 * - Get cracking on all the TODOs spread around
 */

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

  // Private state variables
  AdjacencyList _adjacencies;

  /* Private member functions */
  bool _validAtomIndex(const AtomIndexType& a) const;

public:
/* Public members */
  StereocenterList stereocenters;


/* Constructors */
  // From two elements and a shared bond (min to be considered a Molecule)
  Molecule(
    const Delib::ElementType& a,
    const Delib::ElementType& b,
    const BondType& bondType
  ); 

  // From internal components
  Molecule(const AdjacencyList& adjacencies);
  Molecule(
    const AdjacencyList& adjacencies,
    const StereocenterList& stereocenters
  );

/* Modifiers */
  // Add an atom bonded to an existing atom
  AtomIndexType addAtom(
    const Delib::ElementType& elementType,
    const AtomIndexType& bondedToIndex,
    const BondType& bondType
  );

  // Add a bond betwen existing atoms in the Molecule
  void addBond(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const BondType& bondType
  ); 

  void changeElementType(
    const AtomIndexType& a,
    const Delib::ElementType& elementType
  );

  // Remove an atom. This removes all bonds to and from this atom index.
  void removeAtom(const AtomIndexType& a); 

  // Remove a bond between two atoms.
  void removeBond(
    const AtomIndexType& a,
    const AtomIndexType& b
  );

  void updateStereocenters();

/* Information retrieval */
  void dumpGraphviz(const std::string& filename) const;

  int formalCharge(const AtomIndexType& a) const; // TODO not implemented

  const AdjacencyList& getAdjacencyList() const; 

  // Returns a vector of persistent external indices that the atom is bonded to.
  std::vector<AtomIndexType> getBondedAtomIndices(
    const AtomIndexType& a
  ) const;

  // Get the bond type for specified external indices
  BondType getBondType(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;

  std::vector<AdjacencyList::ExplicitEdge> getEdges() const;

  Delib::ElementType getElementType(
    const AtomIndexType& a
  ) const;

  unsigned getNumAtoms() const;
  unsigned getNumBonds() const;


  unsigned hydrogenCount(const AtomIndexType& a) const;

  int oxidationState(const AtomIndexType& a) const; // TODO not implemented

  std::pair<
    std::vector<AtomIndexType>, // the sorted list of substituent priorities
    std::set< // a set of pairs of AtomIndexTypes that have equal priority
      std::pair<
        AtomIndexType,
        AtomIndexType
      >
    >
  > rankPriority(
    const AtomIndexType& a,
    const std::vector<AtomIndexType>& excludeAdjacent = {}
  ) const;

  /* Testing */
  std::pair<bool, std::string> validate() const noexcept; // TODO not implemented

/* Operators */
  /* An efficient implementation of the following two is imperative.
   * Some ideas for fast differentiation can probably be found from the 
   * wikipedia category "Graph invariants"
   */
  bool operator == (const Molecule& b) const; // TODO not implemented
  bool operator != (const Molecule& b) const; // TODO not implemented

/* Friends */
  /* Output stream operator for easier debugging. */
  friend std::ostream& operator<<(std::ostream& os, const Molecule& mol);
};

}

#endif
