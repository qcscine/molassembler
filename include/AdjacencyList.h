#ifndef INCLUDE_ADJACENCY_LIST_H
#define INCLUDE_ADJACENCY_LIST_H

#include <vector>
#include <algorithm>
#include <cassert>

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
    if(!validate()) {
      throw std::runtime_error(
        "Constructing AdjacencyList from Edges yielded invalid state!"
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
   * Returns the size of the AdjacencyList
   * \returns The size of the AdjacencyList
   */
  AtomIndexType size() const noexcept {
    return _adjacencies.size();
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
    return _adjacencies.at(a);
  }
};

} // eo namespace MoleculeManip

#endif
