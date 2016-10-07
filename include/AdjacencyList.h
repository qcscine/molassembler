#include <vector>
#include <algorithm>
#include <cassert>

#include "common_typedefs.h"

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

public:
/* Public member functions */
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
    _adjacencies[a].emplace_back(b);
    _adjacencies[b].emplace_back(a);
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

/* Operators */
  const std::vector<AtomIndexType>& operator[](const AtomIndexType& a) const {
    assert(_isValidIndex(a));
    return _adjacencies[a];
  }
};

} // eo namespace MoleculeManip
