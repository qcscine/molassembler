#ifndef INCLUDE_EDGES_H
#define INCLUDE_EDGES_H

#include <algorithm>
#include <functional>
#include <boost/optional.hpp>
#include <map>

#include "common_typedefs.h"

/*! @file
 *
 * Contains an intermediate graph representation class to help dealing with
 * small graphs for simple testing.
 */

/* Enhancements:
 * - Could provide custom iterator with .i, .j, .bondType members for easier
 *   for loop iteration. Otherwise, it's .first.first, .first.second, .second,
 *   which is not-so-nice.
 */

namespace MoleculeManip {

/*!
 * A small-time graph representation class whose underlying data is essentially
 * a map between a pair of indices and a bond type. Largely unused except for
 * the construction of small graphs during testing.
 */
class Edges {
public:
  using MapType = std::map<
    std::pair<AtomIndexType, AtomIndexType>,
    BondType
  >;

private:
  MapType _edges;

  std::pair<AtomIndexType, AtomIndexType> _makeOrderedPair(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const {
    return std::pair<AtomIndexType, AtomIndexType>(
      std::min(a, b),
      std::max(a, b)
    );
  }

public:
  Edges() = default;
  Edges(
    const std::initializer_list<
      std::pair<
        MapType::key_type,
        MapType::mapped_type
      >
    >& initializerList
  ) {
    for(const auto& element : initializerList) {
      add(element.first.first, element.first.second, element.second);
    }
  }

  void add(
    const AtomIndexType& a,
    const AtomIndexType& b,
    const BondType& bty
  ) noexcept {
    if(_edges.count(_makeOrderedPair(a, b)) == 0) {
      _edges[_makeOrderedPair(a, b)] = bty;
    }
  }

  void clear() noexcept {
    _edges.clear();
  }

  void remove(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) {
    auto it = _edges.find(_makeOrderedPair(a, b));
    if(it != _edges.end()) {
      _edges.erase(it);
    } else { // fail on error!
      throw std::logic_error("No such edge in list!");
    }
  }

  boost::optional<BondType> get(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const {
    if(_edges.count(_makeOrderedPair(a, b)) == 1) {
      return _edges.at(_makeOrderedPair(a, b));
    } else return boost::none;
  }

  /* Decrements all indices in the map if they are bigger than the specified 
   * invalidated index. This helps preserve contiguity of internal indices.
   * 
   * The implementation is O(N) with bad constants. Creates another map,
   * transforms the keys according to the decrement and overwrites the internal
   * map with the new one.  C++17 would be very helpful with extract, allowing
   * individual modification of nodes that have indices where the indices are
   * larger than the removed atom. That ought to be faster.
   */
  void indexInvalidationUpdate(const AtomIndexType& invalidated) {
    // Re-create edges one-by-one (by copy)
    MapType updated;

    auto decrementIfLarger = [&invalidated](const AtomIndexType& i) -> AtomIndexType {
      if(i > invalidated) return i - 1;
      else return i;
    };

    for(const auto& setIterPair : _edges) {
      updated[
        _makeOrderedPair(
          decrementIfLarger(setIterPair.first.first),
          decrementIfLarger(setIterPair.first.second)
        )
      ] = setIterPair.second;
    }

    // overwrite with updated edge list
    _edges = std::move(updated);
  }

  unsigned size() const noexcept {
    return _edges.size();
  }

  MapType::const_iterator begin() const {
    return _edges.cbegin();
  }

  MapType::const_iterator end() const {
    return _edges.cend();
  }

};

} // namespace MoleculeManip

#endif
