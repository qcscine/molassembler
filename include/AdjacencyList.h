#ifndef INCLUDE_ADJACENCY_LIST_H
#define INCLUDE_ADJACENCY_LIST_H

#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>

#include "Edges.h"

namespace MoleculeManip {

class AdjacencyList {
private:
  /* Private members */
  std::vector<
    std::vector<
      AtomIndexType
    >
  > _adjacencies;

public:
/* Public member functions */
  /* Constructors */
  AdjacencyList() = default;
  AdjacencyList(
    const Edges& edges
  ) {
    for(const auto& edge: edges) {
      // resize if indices do not fit
      if(edge.first.second >= _adjacencies.size()) {
        _adjacencies.resize(edge.first.second + 1);
      }
      addAdjacency(edge.first.first, edge.first.second);
    }
  }

  /* Modification */
  /*!
   * Adds an empty vector for a new vertex to the underlying data 
   * structure.
   */
  AtomIndexType addSlot() noexcept {
    _adjacencies.emplace_back();
    return _adjacencies.size() - 1;
  }

  /*!
   * Adds the passed indices to each other's adjacency lists.
   * \param a The first index
   * \param b The second index
   */
  void addAdjacency(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) noexcept {
    _adjacencies.at(a).emplace_back(b);
    _adjacencies.at(b).emplace_back(a);
  }

  /*!
   * Clears the entire list.
   */
  void clear() noexcept {
    _adjacencies.clear();
  }

  /*!
   * Update the AdjacencyList when an atom index is invalidated to return to 
   * contiguous internal indices.
   */
  void indexInvalidationUpdate(const AtomIndexType& invalidatedIndex) {
    // The according row in the AdjacencyList should already be empty, check
    assert(_adjacencies.at(invalidatedIndex).size() == 0);

    // erase it
    _adjacencies.erase(_adjacencies.begin() + invalidatedIndex);

    // update all indices in the AdjacencyList
    for(auto& adjacencyRow : _adjacencies) {
      std::transform(
        adjacencyRow.begin(),
        adjacencyRow.end(),
        adjacencyRow.begin(), // in-place transform
        [&invalidatedIndex] (const auto& index) {
          if(index > invalidatedIndex) return index - 1;
          else return index;
        }
      );
    }
  }

  /*!
   * Removes an adjacency.
   * \param a The first index
   * \param b The second index
   */
  void removeAdjacency(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) {
    _adjacencies.at(a).erase(
      std::remove(
        _adjacencies.at(a).begin(),
        _adjacencies.at(a).end(),
        b
      )
    );

    _adjacencies.at(b).erase(
      std::remove(
        _adjacencies.at(b).begin(),
        _adjacencies.at(b).end(),
        a
      )
    );
  }

  /* Information */
  //! Dump graphviz source file text representing the AdjacencyList
  std::string dumpGraphviz() const {
    std::stringstream ss;
    ss << "graph G {\n  graph [fontname = \"Arial\", layout = neato];\n"
      << "  node [fontname = \"Arial\", shape = circle, style = filled, fontsize=10];\n"
      << "  edge [fontname = \"Arial\"];\n";

    // Node names
    ss << "node [fillcolor=white, fontcolor=black, width=.3, fixedsize=true]\n";
    for(unsigned i = 0; i < _adjacencies.size(); i++) {
      ss << "\"" << i << "\"";
      if(i != _adjacencies.size() - 1) {
        ss << " ";
      };
    }
    ss << ";\n\n";

    // All forward edges (i.e. first < second)
    for(unsigned i = 0; i < _adjacencies.size(); i++) {
      for(unsigned j = 0; j < _adjacencies[i].size(); j++) {
        if(i < _adjacencies[i][j]) {
          ss << "\"" << i << "\" -- \"" << _adjacencies[i][j] << "\"\n";
        }
      }
    }

    ss << "}";

    return ss.str();
  }

  /*!
   * Get a list of adjacencies for a specified index
   * \param a The index to get
   * \returns The list of adjacencies for that index.
   */
  std::vector<AtomIndexType> getAdjacencies(
    const AtomIndexType& a
  ) const {
    return _adjacencies.at(a);
  }

  /*!
   * Checks whether two vertices are adjacent.
   * \param a The first index
   * \param b The second index
   * \returns whether the two vertices are adjacent.
   */
  bool isAdjacent(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const noexcept {
    return std::find(
      _adjacencies.at(a).begin(),
      _adjacencies.at(a).end(),
      b
    ) != _adjacencies.at(a).end();
  }

  /*!
   * Resize to N
   */
  void resize(const unsigned& N) {
    _adjacencies.resize(N, std::vector<AtomIndexType>());
  }

  /*!
   * Returns the size of the AdjacencyList
   * \returns The size of the AdjacencyList
   */
  AtomIndexType size() const noexcept {
    return _adjacencies.size();
  }

/* Operators */
  const std::vector<AtomIndexType>& operator[](const AtomIndexType& a) const {
    return _adjacencies.at(a);
  }

/* Friends */
  friend struct AdjacencyListValidator;
};

} // eo namespace MoleculeManip

#endif
