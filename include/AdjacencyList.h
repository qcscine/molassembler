#ifndef INCLUDE_ADJACENCY_LIST_H
#define INCLUDE_ADJACENCY_LIST_H

#include <vector>
#include <algorithm>
#include <cassert>

#include "common_typedefs.h"
#include "EdgeList.h"

namespace MoleculeManip {

class AdjacencyList {
private:
  /* Private members */
  std::vector<
    std::vector<
      AtomIndexType
    >
  > _adjacencies;

  /* Private member functions */
  /*!
   * Ensures the passed index is in range.
   * \param a The index to check
   * \returns Whether the index is in range.
   */
  bool _isValidIndex(
    const AtomIndexType& a
  ) const {
    return a < _adjacencies.size();
  }

  bool _areValidIndices(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const {
    return (
      _isValidIndex(a)
      && _isValidIndex(b)
      && a != b
    );
  }

public:
/* Public member functions */
  /* Constructors */
  AdjacencyList() = default;
  AdjacencyList(
    const EdgeList& edges
  ) {
    for(const auto& edge: edges) {
      // resize if indices do not fit
      if(std::max(edge.i, edge.j) >= _adjacencies.size()) {
        _adjacencies.resize(std::max(edge.i, edge.j) + 1);
      }
      addAdjacency(edge.i, edge.j);
    }
    if(!validate()) {
      throw std::runtime_error(
        "Constructing AdjacencyList from EdgeList yielded invalid state!"
      );
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
    assert(_areValidIndices(a, b));

    _adjacencies[a].emplace_back(b);
    _adjacencies[b].emplace_back(a);
  }

  /*!
   * Clears the entire list.
   */
  void clear() noexcept {
    _adjacencies.clear();
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
    assert(_areValidIndices(a, b));

    _adjacencies[a].erase(
      std::remove(
        _adjacencies[a].begin(),
        _adjacencies[a].end(),
        b
      )
    );
    _adjacencies[b].erase(
      std::remove(
        _adjacencies[b].begin(),
        _adjacencies[b].end(),
        a
      )
    );
  }

  /* Information */
  /*!
   * Returns the size of the AdjacencyList
   * \returns The size of the AdjacencyList
   */
  AtomIndexType size() const noexcept {
    return _adjacencies.size();
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
    assert(_areValidIndices(a, b));

    return std::find(
      _adjacencies[a].begin(),
      _adjacencies[a].end(),
      b
    ) != _adjacencies[a].end();
  }

  /*!
   * Get a list of adjacencies for a specified index
   * \param a The index to get
   * \returns The list of adjacencies for that index.
   */
  std::vector<AtomIndexType> getAdjacencies(
    const AtomIndexType& a
  ) const {
    assert(_isValidIndex(a));
    return _adjacencies[a];
  }

  /*!
   * Checks if the current state of the list is valid
   */
  bool validate() const {
    // There is an inverse for every forward bond type
    for(unsigned i = 0; i < _adjacencies.size(); i++) {
      auto adjacents = getAdjacencies(i);
      for(const auto& adjacency: adjacents) {
        if(adjacency >= _adjacencies.size()) return false;
        if(!isAdjacent(adjacency, i)) {
          return false;
        }
      }
    }

    return true;
  }

/* Operators */
  const std::vector<AtomIndexType>& operator[](const AtomIndexType& a) const {
    assert(_isValidIndex(a));
    return _adjacencies[a];
  }
};

} // eo namespace MoleculeManip

#endif
