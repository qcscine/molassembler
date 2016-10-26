#ifndef INCLUDE_EDGE_LIST_H
#define INCLUDE_EDGE_LIST_H

#include <algorithm>
#include <vector>
#include <functional>
#include <cassert>
#include <experimental/optional>

#include "Edge.h"

namespace MoleculeManip {

/*!
 * An ordered List with special lookup functions for Edges.
 */
class EdgeList {
private:
/* Private members */
  std::vector<Edge> _edges;

  /*!
   * Binary searches an ordered vector of edge tuples for a specified 
   * edge, returning whether it was found and the position in the 
   * container it was found at.
   * \param[in] a The first atom index of the edge
   * \param[in] b The second atom index of the edge
   * \returns A pair containing whether it was found, and the position in
   *  the container it was found at. 
   * \pre a, b must fulfill a < b
   * \note Is O( log(N) ), where
   * - N is the number of edges stored in the container
   */
  std::experimental::optional<EdgeIndexType> _binarySearch(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const noexcept {
    /* important safeguard, otherwise
     *  unsigned_type right = _edges.size() - 1
     * is an underflow
     */
    if(_edges.size() == 0) return {};

    // initialize L, M, R
    unsigned_type left = 0;
    unsigned_type right = _edges.size() - 1;
    unsigned_type middle;

    while(left <= right) {
      middle = std::floor(
        ((double) left + (double) right) / 2.0
      );

      if(_edges[middle].smallerThan(a, b)) {
        left = middle + 1;
        continue;
      } else if(_edges[middle].greaterThan(a, b)) {
        // underflow safeguard
        if(middle == 0) right = middle;
        else right = middle - 1;
        continue;
      } else { // not smaller and not greater -> equals
        // found it!
        return middle;
      }
    }

    // search fails
    return {};
  }

  /*!
   * Binary searches the ordered vector of edge tuples and returns an
   * index before which the specified edge can be inserted.
   * \param[in] a The first atom index of the edge
   * \param[in] b The second atom index of the edge
   * \returns The position within the container before which the
   *  specified edge can be inserted.
   * \pre a, b must fulfill a < b
   * \note Is O( log(N) ), where
   * - N is the number of edges stored in the container
   */
  EdgeIndexType _binarySearchInsertPosition(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const noexcept {
    /* important safeguard, otherwise
     *  unsigned_type right = container.size() - 1;
     * causes an underflow
     */
    if(_edges.size() == 0) {
      return 0;
    }

    // initializations
    unsigned_type left = 0;
    unsigned_type right = _edges.size() - 1;
    unsigned_type middle;

    while(left <= right) {

      middle = std::floor(
        ((double) left + (double) right) / 2.0
      );

      if(left == right) {
        if(_edges[middle].smallerThan(a, b)) {
          return middle + 1;
        } else {
          return middle;
        }
      }

      if(_edges[middle].smallerThan(a, b)) {
        left = middle + 1;
        continue;
      }

      if(_edges[middle].greaterThan(a, b)) {
        // underflow safeguard
        if(middle == 0) right = middle;
        else right = middle - 1; 
        continue;
      }
    }

    // "majority vote"
    if(middle == left) {
      return left;
    } else if(middle == right) {
      return right;
    } else {
      return 0;
    }
  }

  void _remove(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) {
    auto foundOption = _binarySearch(
      a,
      b
    );
    if(foundOption) {
      _edges.erase(
        _edges.begin() + foundOption.value()
      );
    }
  }

public:
/* Constructors */
  EdgeList() = default;
  EdgeList(const std::vector<Edge>& edges) {
    for(const auto& edge : edges) {
      add(edge);
    }
  }

/* Modification */
  void add(const Edge& edge) noexcept {
    auto insert_position = _binarySearchInsertPosition(
      edge.i,
      edge.j
    );

    _edges.insert(
      _edges.begin() + insert_position,
      edge
    );
  }

  void clear() noexcept {
    _edges.clear();
  }

  void remove(const Edge& edge) noexcept {
    _remove(edge.i, edge.j);
  }

  void remove(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) noexcept {
    _remove(
      std::min(a, b),
      std::max(a, b)
    );
  }

/* Information */
  inline Edge get(const EdgeIndexType& position) const noexcept {
    assert(position < _edges.size());
    return _edges[position];
  }

  /*!
   * Checks whether a list of EdgeTypes is ordered.
   * \param 
   * \returns whether the list is ordered.
   * \note is O(E) 
   */
  bool is_ordered() const noexcept {
    /* this is a necessary underflow guard:
     * in the loop condition all types are unsigned:
     *  i < _edges.size() - 1
     * if _edges.size() is 0, the expression is an underflow
     */
    if(_edges.size() == 0) return true;

    for(EdgeIndexType i = 0; i < _edges.size() - 1; i++) {
      if(!(
        _edges[i] < _edges[i + 1]
      )) {
        return false;
      }
    }

    return true;
  }

  std::experimental::optional<EdgeIndexType> search(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const noexcept {
    return _binarySearch(
      std::min(a, b),
      std::max(a, b)
    );
  }

  unsigned size() const noexcept {
    return _edges.size();
  }

/* Iterators */
  std::vector<Edge>::const_iterator begin() const {
    return _edges.cbegin();
  }
  std::vector<Edge>::const_iterator end() const {
    return _edges.cend();
  }
}; 

} // eo namespace

#endif
