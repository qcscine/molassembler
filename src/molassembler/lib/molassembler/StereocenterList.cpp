#include "molassembler/StereocenterList.h"

#include "temple/Containers.h"

namespace molassembler {

void StereocenterList::add(
  AtomIndexType i,
  AtomStereocenter stereocenter
) {
  _atomStereocenters.emplace(
    i,
    std::move(stereocenter)
  );
}

void StereocenterList::add(
  GraphType::edge_descriptor edge,
  BondStereocenter stereocenter
) {
  _bondStereocenters.emplace(
    edge,
    std::move(stereocenter)
  );
}

void StereocenterList::clear() {
  _atomStereocenters.clear();
  _bondStereocenters.clear();
}

void StereocenterList::clearBonds() {
  _bondStereocenters.clear();
}

void StereocenterList::propagateVertexRemoval(const AtomIndexType removedIndex) {
  // Drop any stereocenters involving this atom from the list
  auto findIter = _atomStereocenters.find(removedIndex);
  if(findIter != _atomStereocenters.end()) {
    _atomStereocenters.erase(findIter);
  }

  /* Go through all state in the StereocenterList and decrement any indices
   * larger than the one being removed
   */
  for(auto& stereocenters : atomStereocenters()) {
    stereocenters.propagateVertexRemoval(removedIndex);
  }

  /*for(auto& bondMapPair : bondStereocenters()) {
    bondMapPair.second.propagateVertexRemoval(removedIndex);
  }*/
}

void StereocenterList::remove(const AtomIndexType index) {
  auto findIter = _atomStereocenters.find(index);

  if(findIter != _atomStereocenters.end()) {
    _atomStereocenters.erase(findIter);
  } else {
    throw std::logic_error("No such atom stereocenter found!");
  }
}

void StereocenterList::remove(const GraphType::edge_descriptor edge) {
  auto findIter = _bondStereocenters.find(edge);

  if(findIter != _bondStereocenters.end()) {
    _bondStereocenters.erase(findIter);
  } else {
    throw std::logic_error("No such bond stereocenter found!");
  }
}

void StereocenterList::try_remove(const AtomIndexType index) {
  auto findIter = _atomStereocenters.find(index);

  if(findIter != _atomStereocenters.end()) {
    _atomStereocenters.erase(findIter);
  }
}

void StereocenterList::try_remove(const GraphType::edge_descriptor edge) {
  auto findIter = _bondStereocenters.find(edge);

  if(findIter != _bondStereocenters.end()) {
    _bondStereocenters.erase(findIter);
  }
}

boost::optional<AtomStereocenter&> StereocenterList::option(const AtomIndexType index) {
  auto findIter = _atomStereocenters.find(index);

  if(findIter != _atomStereocenters.end()) {
    return findIter->second;
  }

  return boost::none;
}

boost::optional<BondStereocenter&> StereocenterList::option(const GraphType::edge_descriptor edge) {
  auto findIter = _bondStereocenters.find(edge);

  if(findIter != _bondStereocenters.end()) {
    return findIter->second;
  }

  return boost::none;
}

/* Information */
bool StereocenterList::empty() const {
  return _atomStereocenters.empty() && _bondStereocenters.empty();
}

unsigned StereocenterList::size() const {
  return _atomStereocenters.size() + _bondStereocenters.size();
}

boost::optional<const AtomStereocenter&> StereocenterList::option(const AtomIndexType index) const {
  auto findIter = _atomStereocenters.find(index);

  if(findIter != _atomStereocenters.end()) {
    return findIter->second;
  }

  return boost::none;
}

boost::optional<const BondStereocenter&> StereocenterList::option(const GraphType::edge_descriptor edge) const {
  auto findIter = _bondStereocenters.find(edge);

  if(findIter != _bondStereocenters.end()) {
    return findIter->second;
  }

  return boost::none;
}


/* Iterators */
boost::range_detail::select_second_mutable_range<StereocenterList::AtomMapType> StereocenterList::atomStereocenters() {
  return _atomStereocenters | boost::adaptors::map_values;
}

boost::range_detail::select_second_const_range<StereocenterList::AtomMapType> StereocenterList::atomStereocenters() const {
  return _atomStereocenters | boost::adaptors::map_values;
}

boost::range_detail::select_second_mutable_range<StereocenterList::BondMapType> StereocenterList::bondStereocenters() {
  return _bondStereocenters | boost::adaptors::map_values;
}

boost::range_detail::select_second_const_range<StereocenterList::BondMapType> StereocenterList::bondStereocenters() const {
  return _bondStereocenters | boost::adaptors::map_values;
}

bool StereocenterList::hasZeroAssignmentStereocenters() const {
  return temple::any_of(
    atomStereocenters(),
    [](const auto& stereocenter) -> bool {
      return stereocenter.numAssignments() == 0u;
    }
  ) || temple::any_of(
    bondStereocenters(),
    [](const auto& stereocenter) -> bool {
      return stereocenter.numAssignments() == 0u;
    }
  );
}

bool StereocenterList::hasUnassignedStereocenters() const {
  return temple::any_of(
    atomStereocenters(),
    [](const auto& stereocenter) -> bool {
      return !stereocenter.assigned();
    }
  ) || temple::any_of(
    bondStereocenters(),
    [](const auto& stereocenter) -> bool {
      return !stereocenter.assigned();
    }
  );
}

bool StereocenterList::compare(
  const StereocenterList& other,
  const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
) const {
  if(comparisonBitmask & AtomEnvironmentComponents::Symmetries) {
    // Check sizes
    if(
      _atomStereocenters.size() != other._atomStereocenters.size()
      || _bondStereocenters.size() != other._bondStereocenters.size()
    ) {
      return false;
    }

    // Check all atom stereocenters
    for(const auto& stereocenter : _atomStereocenters | boost::adaptors::map_values) {
      auto otherStereocenterOption = other.option(stereocenter.centralIndex());

      // Ensure there is a matching stereocenter
      if(!otherStereocenterOption) {
        return false;
      }

      // Ensure the symmetries match
      if(stereocenter.getSymmetry() != otherStereocenterOption->getSymmetry()) {
        return false;
      }

      if(comparisonBitmask & AtomEnvironmentComponents::Stereopermutations) {
        // Ensure stereopermutations match
        if(stereocenter.assigned() != otherStereocenterOption->assigned()) {
          return false;
        }
      }
    }

    // Check all bond stereocenters
    for(const auto& stereocenter : _bondStereocenters | boost::adaptors::map_values) {
      auto otherStereocenterOption = other.option(stereocenter.edge());

      // Ensure there is a matching stereocenter
      if(!otherStereocenterOption) {
        return false;
      }

      // Ensure construction state matches
      if(!stereocenter.hasSameCompositeOrientation(otherStereocenterOption.value())) {
        return false;
      }

      if(comparisonBitmask & AtomEnvironmentComponents::Stereopermutations) {
        // Ensure stereopermutations match
        if(stereocenter.assigned() != otherStereocenterOption->assigned()) {
          return false;
        }
      }
    }
  }

  return true;
}

/* Operators */
bool StereocenterList::operator == (const StereocenterList& other) const {
  return (
    _atomStereocenters == other._atomStereocenters
    && _bondStereocenters == other._bondStereocenters
  );
}

bool StereocenterList::operator != (const StereocenterList& other) const {
  return !(*this == other);
}

} // namespace molassembler
