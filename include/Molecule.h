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
#include "Edges.h"
#include "StereocenterList.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "BondDistance.h"

namespace MoleculeManip {

/* TODO
 * - Get cracking on all the TODOs
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

  // The set of QC data on the atoms
  Delib::ElementTypeCollection _elements;

  // The information on interconnectedness of the atoms
  AdjacencyList _adjacencies;
  Edges _edges;
  StereocenterList _stereocenters;

  /* Private member functions */
  void _detectStereocenters();
  void _dumpGraphviz(std::ostream& os) const;
  std::vector<DistanceConstraint> _createConstraint(
    const std::vector<AtomIndexType>& chain
  ) const;
  bool _validAtomIndex(const AtomIndexType& a) const;

public:
/* Constructors */
  // From two elements and a shared bond (min to be considered a Molecule)
  Molecule(
    const Delib::ElementType& a,
    const Delib::ElementType& b,
    const BondType& bondType
  ); 

  // From the internal components
  Molecule(
    const Delib::ElementTypeCollection& elements,
    const AdjacencyList& adjacencies,
    const Edges& edges
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

  // Remove an atom. This removes all bonds to and from this atom index.
  void removeAtom(const AtomIndexType& a); 

  // Remove a bond between two atoms.
  void removeBond(
    const AtomIndexType& a,
    const AtomIndexType& b
  );

/* Information retrieval */

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

  // Get a vector of Distance Geometry chirality constraints
  std::vector<ChiralityConstraint> getChiralityConstraints() const;

  // Get a matrix containing Distance Geometry distance bounds
  DistanceGeometry::DistanceBoundsMatrix getDistanceBoundsMatrix() const;

  const Edges& getEdges() const; 

  Delib::ElementType getElementType(
    const AtomIndexType& a
  ) const;

  AtomIndexType getNumAtoms() const;
  EdgeIndexType getNumBonds() const;

  unsigned hydrogenCount(const AtomIndexType& a) const;

  int oxidationState(const AtomIndexType& a) const; // TODO not implemented

  std::pair<
    std::vector<AtomIndexType>, // the sorted list of substituent priorities
    std::set< // a set of pairs of AtomIndexTypes that are EQUAL
      std::pair<
        AtomIndexType,
        AtomIndexType
      >
    >
  > rankPriority(
    const AtomIndexType& a,
    const std::vector<AtomIndexType>& excludeAdjacent
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
