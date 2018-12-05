// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/StereopermutatorList.h"

#include "temple/Functional.h"

namespace molassembler {

void StereopermutatorList::add(
  AtomIndex i,
  AtomStereopermutator stereopermutator
) {
  _atomStereopermutators.emplace(
    i,
    std::move(stereopermutator)
  );
}

void StereopermutatorList::add(
  const BondIndex& edge,
  BondStereopermutator stereopermutator
) {
  _bondStereopermutators.emplace(
    edge,
    std::move(stereopermutator)
  );
}

void StereopermutatorList::clear() {
  _atomStereopermutators.clear();
  _bondStereopermutators.clear();
}

void StereopermutatorList::clearBonds() {
  _bondStereopermutators.clear();
}

void StereopermutatorList::propagateVertexRemoval(const AtomIndex removedIndex) {
  // Drop any stereopermutators involving this atom from the list
  auto findIter = _atomStereopermutators.find(removedIndex);
  if(findIter != _atomStereopermutators.end()) {
    _atomStereopermutators.erase(findIter);
  }

  /* Go through all state in the StereopermutatorList and decrement any indices
   * larger than the one being removed
   */
  for(auto& stereopermutators : atomStereopermutators()) {
    stereopermutators.propagateVertexRemoval(removedIndex);
  }

  /*for(auto& bondMapPair : bondStereopermutators()) {
    bondMapPair.second.propagateVertexRemoval(removedIndex);
  }*/
}

void StereopermutatorList::remove(const AtomIndex index) {
  auto findIter = _atomStereopermutators.find(index);

  if(findIter != _atomStereopermutators.end()) {
    _atomStereopermutators.erase(findIter);
  } else {
    throw std::logic_error("No such atom stereopermutator found!");
  }
}

void StereopermutatorList::remove(const BondIndex& edge) {
  auto findIter = _bondStereopermutators.find(edge);

  if(findIter != _bondStereopermutators.end()) {
    _bondStereopermutators.erase(findIter);
  } else {
    throw std::logic_error("No such bond stereopermutator found!");
  }
}

void StereopermutatorList::try_remove(const AtomIndex index) {
  auto findIter = _atomStereopermutators.find(index);

  if(findIter != _atomStereopermutators.end()) {
    _atomStereopermutators.erase(findIter);
  }
}

void StereopermutatorList::try_remove(const BondIndex& edge) {
  auto findIter = _bondStereopermutators.find(edge);

  if(findIter != _bondStereopermutators.end()) {
    _bondStereopermutators.erase(findIter);
  }
}

boost::optional<AtomStereopermutator&> StereopermutatorList::option(const AtomIndex index) {
  auto findIter = _atomStereopermutators.find(index);

  if(findIter != _atomStereopermutators.end()) {
    return findIter->second;
  }

  return boost::none;
}

boost::optional<BondStereopermutator&> StereopermutatorList::option(const BondIndex& edge) {
  auto findIter = _bondStereopermutators.find(edge);

  if(findIter != _bondStereopermutators.end()) {
    return findIter->second;
  }

  return boost::none;
}

/* Information */
bool StereopermutatorList::empty() const {
  return _atomStereopermutators.empty() && _bondStereopermutators.empty();
}

unsigned StereopermutatorList::A() const {
  return _atomStereopermutators.size();
}

unsigned StereopermutatorList::B() const {
  return _bondStereopermutators.size();
}

unsigned StereopermutatorList::size() const {
  return _atomStereopermutators.size() + _bondStereopermutators.size();
}

boost::optional<const AtomStereopermutator&> StereopermutatorList::option(const AtomIndex index) const {
  auto findIter = _atomStereopermutators.find(index);

  if(findIter != _atomStereopermutators.end()) {
    return findIter->second;
  }

  return boost::none;
}

boost::optional<const BondStereopermutator&> StereopermutatorList::option(const BondIndex& edge) const {
  auto findIter = _bondStereopermutators.find(edge);

  if(findIter != _bondStereopermutators.end()) {
    return findIter->second;
  }

  return boost::none;
}


/* Iterators */
boost::range_detail::select_second_mutable_range<StereopermutatorList::AtomMapType> StereopermutatorList::atomStereopermutators() {
  return _atomStereopermutators | boost::adaptors::map_values;
}

boost::range_detail::select_second_const_range<StereopermutatorList::AtomMapType> StereopermutatorList::atomStereopermutators() const {
  return _atomStereopermutators | boost::adaptors::map_values;
}

boost::range_detail::select_second_mutable_range<StereopermutatorList::BondMapType> StereopermutatorList::bondStereopermutators() {
  return _bondStereopermutators | boost::adaptors::map_values;
}

boost::range_detail::select_second_const_range<StereopermutatorList::BondMapType> StereopermutatorList::bondStereopermutators() const {
  return _bondStereopermutators | boost::adaptors::map_values;
}

bool StereopermutatorList::hasZeroAssignmentStereopermutators() const {
  return temple::any_of(
    atomStereopermutators(),
    [](const auto& stereopermutator) -> bool {
      return stereopermutator.numAssignments() == 0u;
    }
  ) || temple::any_of(
    bondStereopermutators(),
    [](const auto& stereopermutator) -> bool {
      return stereopermutator.numAssignments() == 0u;
    }
  );
}

bool StereopermutatorList::hasUnassignedStereopermutators() const {
  return temple::any_of(
    atomStereopermutators(),
    [](const auto& stereopermutator) -> bool {
      return !stereopermutator.assigned();
    }
  ) || temple::any_of(
    bondStereopermutators(),
    [](const auto& stereopermutator) -> bool {
      return !stereopermutator.assigned();
    }
  );
}

bool StereopermutatorList::compare(
  const StereopermutatorList& other,
  const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
) const {
  if(comparisonBitmask & AtomEnvironmentComponents::Symmetries) {
    // Check sizes
    if(
      _atomStereopermutators.size() != other._atomStereopermutators.size()
      || _bondStereopermutators.size() != other._bondStereopermutators.size()
    ) {
      return false;
    }

    // Check all atom stereopermutators
    for(const auto& stereopermutator : _atomStereopermutators | boost::adaptors::map_values) {
      auto otherStereopermutatorOption = other.option(stereopermutator.centralIndex());

      // Ensure there is a matching stereopermutator
      if(!otherStereopermutatorOption) {
        return false;
      }

      // Ensure the symmetries match
      if(stereopermutator.getSymmetry() != otherStereopermutatorOption->getSymmetry()) {
        return false;
      }

      if(comparisonBitmask & AtomEnvironmentComponents::Stereopermutations) {
        // Ensure stereopermutations match
        if(stereopermutator.assigned() != otherStereopermutatorOption->assigned()) {
          return false;
        }
      }
    }

    // Check all bond stereopermutators
    for(const auto& stereopermutator : _bondStereopermutators | boost::adaptors::map_values) {
      auto otherStereopermutatorOption = other.option(stereopermutator.edge());

      // Ensure there is a matching stereopermutator
      if(!otherStereopermutatorOption) {
        return false;
      }

      // Ensure construction state matches
      if(!stereopermutator.hasSameCompositeOrientation(otherStereopermutatorOption.value())) {
        return false;
      }

      if(comparisonBitmask & AtomEnvironmentComponents::Stereopermutations) {
        // Ensure stereopermutations match
        if(stereopermutator.assigned() != otherStereopermutatorOption->assigned()) {
          return false;
        }
      }
    }
  }

  return true;
}

/* Operators */
bool StereopermutatorList::operator == (const StereopermutatorList& other) const {
  return (
    _atomStereopermutators == other._atomStereopermutators
    && _bondStereopermutators == other._bondStereopermutators
  );
}

bool StereopermutatorList::operator != (const StereopermutatorList& other) const {
  return !(*this == other);
}

} // namespace molassembler
