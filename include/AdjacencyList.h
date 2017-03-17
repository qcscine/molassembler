#ifndef INCLUDE_ADJACENCY_LIST_H
#define INCLUDE_ADJACENCY_LIST_H

// STL
#include <fstream>

// Libraries
#include "boost/graph/graphviz.hpp"

// Delib
#include "ElementInfo.h"
#include "Types/ElementTypeCollection.h"

#include "common_typedefs.h"
#include "Edges.h"
#include "RangeForTemporary.h"

namespace MoleculeManip {

class AdjacencyList {
private:
/* State */
  GraphType _adjacencies;

/* Members */
  // Helper class to write the Graph as Graphviz output
  struct MolGraphWriter {
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
     {BondType::Single, "color = \"black\""},
      {BondType::Double, "color = \"black:invis:black\""},
      {BondType::Triple, "color = \"black:invis:black:invis:black\""},
      {BondType::Quadruple, "label = \"4\""},
      {BondType::Quintuple, "label = \"5\""},
      {BondType::Sextuple, "label = \"6\""},
      {BondType::Aromatic, "style = \"dashed\""},
      {BondType::Eta, "style = \"dotted\""}
    };

    /* State */
    // We promise to be good and not change anything
    const GraphType* const graphPtr;

    /* Constructor */
    MolGraphWriter(const GraphType* passGraphPtr) : graphPtr(passGraphPtr) {}

    /* Information */
    Delib::ElementType getElementType(const AtomIndexType& vertexIndex) const {
      return (*graphPtr)[vertexIndex].elementType;
    }

    // Global options
    void operator() (std::ostream& os) const {
      os << "graph [fontname = \"Arial\", layout = neato];\n"
        << "node [fontname = \"Arial\", shape = circle, style = filled];\n"
        << "edge [fontname = \"Arial\"];\n";
    }

    // Vertex options
    void operator() (std::ostream& os, const AtomIndexType& vertexIndex) const {
      const std::string symbolString = Delib::ElementInfo::symbol(
        getElementType(vertexIndex)
      );

      os << "[";
      
      // Add element name and index label
      os << "label = \"" << symbolString << vertexIndex << "\"";

      // Coloring
      if(elementBGColorMap.count(symbolString)) {
        os << ", fillcolor=\"" << elementBGColorMap.at(symbolString) << "\"";
      } else { // default
        os << ", fillcolor=\"white\"";
      }
      if(elementTextColorMap.count(symbolString)) {
        os << ", fontcolor=\"" << elementTextColorMap.at(symbolString) << "\"";
      } else { // default
        os << ", fontcolor=\"orange\"";
      }

      // Font sizing
      if(symbolString == "H") os << ", fontsize=10, width=.3, fixedsize=true";
      
      os << "]";
    }

    // Edge options
    void operator() (std::ostream& os, const EdgeIndexType& edgeIndex) const {
      os << "[";

      // Bond Type display options
      auto bondType = (*graphPtr)[edgeIndex].bondType;
      if(bondTypeDisplayString.count(bondType)) {
        os << bondTypeDisplayString.at(bondType);
      }

      // If one of the bonded atoms is a hydrogen, shorten the bond
      if(
        getElementType(
          boost::target(edgeIndex, *graphPtr)
        ) == Delib::ElementType::H
        || getElementType(
          boost::source(edgeIndex, *graphPtr)
        ) == Delib::ElementType::H
      ) {
        os << ", len=0.5";
      }

      os << "]";
    }
  };

/* Private members */
  bool _isValidIndex(const AtomIndexType& index) const {
    return index < size();
  }

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
  ) {
    for(const auto& element: elements) {
      addAtom(element);
    }

    for(const auto& edge: edges) {
      addBond(edge.first.first, edge.first.second, edge.second);
    }
  }

/* Modification */
  AtomIndexType addAtom(
    const Delib::ElementType& elementType
  ) {
    auto vertex = boost::add_vertex(_adjacencies);
    _adjacencies[vertex].elementType = elementType;
    return vertex;
  }

  void addBond(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const BondType& bondType
  ) {
    assert(_isValidIndex(a) && _isValidIndex(b) && a != b);
    auto edgeAddPair = boost::add_edge(a, b, _adjacencies);
    _adjacencies[edgeAddPair.first].bondType = bondType;
  }

  void clear() {
    // Delete EVERYTHING
    _adjacencies.clear();
  }

  void removeAtom(const AtomIndexType& a) {
    // Remove all edges to and from this vertex
    boost::clear_vertex(a, _adjacencies);

    // Remove the vertex itself
    boost::remove_vertex(a, _adjacencies);
  }

  void removeBond(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) {
    // Find edges
    auto edgePair = boost::edge(a, b, _adjacencies);
    if(edgePair.second) {
      boost::remove_edge(edgePair.first, _adjacencies);
    }
  }

/* Information */
  const GraphType& access() const {
    return _adjacencies;
  }

  bool isAdjacent(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const {
    GraphType::adjacency_iterator begin, end;
    std::tie(begin, end) = boost::adjacent_vertices(a, _adjacencies);
    return std::find(
      begin,
      end,
      b
    ) != end;
  }

  unsigned getNumAdjacencies(
    const AtomIndexType& a
  ) const {
    return boost::out_degree(a, _adjacencies);
  }

  unsigned getNumNonEtaAdjacencies(
    const AtomIndexType& a
  ) const {
    unsigned count = 0;

    for(
      const auto& edgeIndex:
      RangeForTemporary<GraphType::out_edge_iterator>(
        boost::out_edges(a, _adjacencies)
      )
    ) {
      if(_adjacencies[edgeIndex].bondType != BondType::Eta) {
        count += 1;
      }
    }

    return count;
  }

  /*! Returns a range-for temporary object allowing c++11 style for loop 
   * iteration through an atom's adjacencies
   */
  RangeForTemporary<GraphType::adjacency_iterator> iterateAdjacencies(
    const AtomIndexType& a
  ) const {
    return RangeForTemporary<
      GraphType::adjacency_iterator
    >(
      boost::adjacent_vertices(a, _adjacencies)
    );
  }

  // Creates a copy of the contained data suitable for the Edges class
  std::vector<ExplicitEdge> getEdges() const {
    std::vector<ExplicitEdge> edges;

    for(
      const auto& edgeIndex : 
      RangeForTemporary<
        GraphType::edge_iterator
      >(boost::edges(_adjacencies))
    ) {
      edges.push_back(ExplicitEdge({
        {boost::source(edgeIndex, _adjacencies), boost::target(edgeIndex, _adjacencies)},
        _adjacencies[edgeIndex].bondType
      }));
    }

    return edges;
  }

  std::vector<AtomIndexType> getAdjacencies(
    const AtomIndexType& a
  ) const {
    std::vector<AtomIndexType> copy;

    // C++17 auto [begin, end] = ...
    GraphType::adjacency_iterator begin, end;
    std::tie(begin, end) = boost::adjacent_vertices(a, _adjacencies);
    std::copy(
      begin,
      end,
      std::back_inserter(copy)
    );

    return copy;
  }

  Delib::ElementType getElementType(const AtomIndexType& index) const {
    assert(_isValidIndex(index));
    return _adjacencies[index].elementType;
  }

  boost::optional<BondType> getBondType(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const {
    auto edgePair = boost::edge(a, b, _adjacencies);

    if(edgePair.second) {
      return _adjacencies[edgePair.first].bondType;
    } else return boost::none;

  }

  // TODO deprecate and remove
  unsigned size() const noexcept {
    return boost::num_vertices(_adjacencies);
  }

  unsigned nAtoms() const {
    return boost::num_vertices(_adjacencies);
  }

  unsigned nBonds() const {
    return boost::num_edges(_adjacencies);
  }

  void dumpGraphviz(const std::string& filename) const {
    MolGraphWriter propertyWriter(&_adjacencies);

    std::ofstream outStream(filename);

    boost::write_graphviz(
      outStream,
      _adjacencies,
      propertyWriter,
      propertyWriter,
      propertyWriter
    );

    outStream.close();
  }

/* Operators */
  RangeForTemporary<GraphType::adjacency_iterator> operator[](
    const AtomIndexType& a
  ) const {
    return RangeForTemporary<
      GraphType::adjacency_iterator
    >(
      boost::adjacent_vertices(a, _adjacencies)
    );
  }

/* Friends */
  friend struct AdjacencyListValidator;
};

} // eo namespace MoleculeManip

#endif
