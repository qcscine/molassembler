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
#include "VSEPR.h"

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

  // State variables
  Delib::ElementTypeCollection _elements;
  AdjacencyList _adjacencies;
  Edges _edges;
  StereocenterList _stereocenters;

  /* Private member functions */
  std::vector<DistanceConstraint> _createConstraint(
    const std::vector<AtomIndexType>& chain
  ) const;
  void _detectStereocenters();
  void _dumpGraphviz(std::ostream& os) const;
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

/* Temporary experimentation */
  std::vector<LocalGeometry::LigandType> _reduceToLigandTypes(
    const AtomIndexType& index
  ) {
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
    assert(_adjacencies.getAdjacencies(index).size() > 1);

    // first basic stuff for VSEPR, later L and X for transition metals
    // geometry inference does not care if the substituents are somehow 
    // connected (unless in later models the entire structure is considered)
    std::vector<LocalGeometry::LigandType> ligands;

    for(const auto& adjacentIndex: _adjacencies.getAdjacencies(index)) {
      ligands.push_back(
        LocalGeometry::LigandType {
          0, 0, {{  // L and X are 0 since only VSEPR is considered for now
            _elements[adjacentIndex],
            getBondType(index, adjacentIndex)
          }}
        }
      );
    }

    return ligands;
  }

  std::map<AtomIndexType, Symmetry::Name> _determineLocalGeometries() {
    std::map<AtomIndexType, Symmetry::Name> symmetryMap;

    for(AtomIndexType i = 0; i < _adjacencies.size(); i++) {
      if(_adjacencies.getAdjacencies(i).size() > 1) {
        auto ligandsVector = _reduceToLigandTypes(i);
        // TODO this below is invalid for metals!
        unsigned nSites = _adjacencies.getAdjacencies(i).size();
        int formalCharge = 0;

        symmetryMap[i] = LocalGeometry::VSEPR::determineGeometry(
          _elements[i],
          nSites,
          ligandsVector,
          formalCharge
        );
      }
    }

    return symmetryMap;
  }

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
