#ifndef INCLUDE_EDGES_H
#define INCLUDE_EDGES_H

#include <algorithm>
#include <functional>
#include <boost/optional.hpp>
#include <map>

#include "common_typedefs.h"

/* Enhancements:
 * - Could provide custom iterator with .i, .j, .bondType members for easier
 *   for loop iteration. Otherwise, it's .first.first, .first.second, .second,
 *   which is not-so-nice.
 */

namespace MoleculeManip {

class Edges {
public:
  using MapType = std::map<
    std::pair<AtomIndexType, AtomIndexType>,
    BondType
  >;
private:
  MapType _edges;

  std::pair<AtomIndexType, AtomIndexType> _makePair(
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
    if(_edges.count(_makePair(a, b)) == 0) {
      _edges[_makePair(a, b)] = bty;
    }
  }

  void clear() noexcept {
    _edges.clear();
  }

  void remove(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) noexcept {
    auto it = _edges.find(_makePair(a, b));
    if(it != _edges.end()) {
      _edges.erase(it);
    }
  }

  boost::optional<BondType> get(
    const AtomIndexType& a,
    const AtomIndexType& b
  ) const {
    if(_edges.count(_makePair(a, b)) == 1) {
      return _edges.at(_makePair(a, b));
    } else return boost::none;
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

} // eo namespace

#endif
