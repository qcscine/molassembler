#ifndef INCLUDE_ADJACENCY_LIST_H
#define INCLUDE_ADJACENCY_LIST_H

// STL
#include <fstream>

// Libraries
#include "boost/graph/graphviz.hpp"

// Delib
#include "Delib/ElementInfo.h"
#include "Delib/ElementTypeCollection.h"
#include "Delib/PositionCollection.h"

#include "common_typedefs.h"
#include "Edges.h"
#include "RangeForTemporary.h"
#include "StereocenterList.h"

#include "symmetry_information/Symmetries.h"
#include "VSEPR.h"

namespace MoleculeManip {

class AdjacencyList {
private:
/* State */
  GraphType _adjacencies;

/* Members */
  // Helper class to write the Graph as Graphviz output
  struct MolGraphWriter;

/* Private members */
  bool _isValidIndex(const AtomIndexType& index) const;
  std::vector<AtomIndexType> _getCNStereocenterCandidates() const;
  std::vector<EdgeIndexType> _getEZStereocenterCandidates() const;
  std::vector<LocalGeometry::LigandType> _reduceToLigandTypes(
    const AtomIndexType& index
  ) const;

public:
/* Typedefs */
  using ExplicitEdge = std::pair<
    Edges::MapType::key_type,
    Edges::MapType::mapped_type
  >;

/* Constructors */
  AdjacencyList() = default;
  AdjacencyList(
    const Delib::ElementTypeCollection& elements,
    const Edges& edges
  );

/* Modification */
  AtomIndexType addAtom(const Delib::ElementType& elementType);

  void addBond(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const BondType& bondType
  );

  void changeElementType(
    const AtomIndexType& a,
    const Delib::ElementType& elementType
  );

  void clear();

  void removeAtom(const AtomIndexType& a);

  void removeBond(
    const AtomIndexType& a,
    const AtomIndexType& b
  );

/* Information */
  const GraphType& access() const;

  bool isAdjacent(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;

  StereocenterList detectStereocenters() const;

  Symmetry::Name determineLocalGeometry(
    const AtomIndexType& index
  ) const;

  unsigned getNumAdjacencies(const AtomIndexType& a) const;

  unsigned getNumNonEtaAdjacencies(const AtomIndexType& a) const;

  // Creates a copy of the contained data suitable for the Edges class
  std::vector<ExplicitEdge> getEdges() const;

  std::vector<AtomIndexType> getAdjacencies(const AtomIndexType& a) const;

  Delib::ElementType getElementType(const AtomIndexType& index) const;

  boost::optional<BondType> getBondType(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const;

  StereocenterList inferStereocentersFromPositions(
    const Delib::PositionCollection& positions
  ) const;

  /*! Returns a range-for temporary object allowing c++11 style for loop 
   * iteration through an atom's adjacencies
   */
  RangeForTemporary<GraphType::adjacency_iterator> iterateAdjacencies(
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
  > rankPriority(
    const AtomIndexType& a,
    const std::vector<AtomIndexType>& excludeAdjacent = {}
  ) const;

  unsigned numAtoms() const;

  unsigned numBonds() const;

  void dumpGraphviz(const std::string& filename) const;

/* Operators */
  RangeForTemporary<GraphType::adjacency_iterator> operator[](
    const AtomIndexType& a
  ) const;

/* Friends */
  friend struct AdjacencyListValidator;
};

} // eo namespace MoleculeManip

#endif
