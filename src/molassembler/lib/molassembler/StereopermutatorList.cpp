/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/StereopermutatorList.h"

#include "temple/Functional.h"

namespace Scine {

namespace molassembler {

StereopermutatorList::AtomMapType::iterator StereopermutatorList::add(
  AtomStereopermutator stereopermutator
) {
  AtomIndex i = stereopermutator.centralIndex();

  auto emplaceResult = _atomStereopermutators.emplace(
    i,
    std::move(stereopermutator)
  );

  if(!emplaceResult.second) {
    throw std::logic_error("Stereopermutator not added. Another is already at its place");
  }

  return emplaceResult.first;
}

StereopermutatorList::BondMapType::iterator StereopermutatorList::add(
  BondStereopermutator stereopermutator
) {
  BondIndex edge = stereopermutator.edge();

  auto emplaceResult = _bondStereopermutators.emplace(
    edge,
    std::move(stereopermutator)
  );

  if(!emplaceResult.second) {
    throw std::logic_error("Stereopermutator not added. Another is already at its place");
  }

  return emplaceResult.first;
}


//! Apply an index mapping to the list of stereopermutators
void StereopermutatorList::applyPermutation(const std::vector<AtomIndex>& permutation) {
  // C++17 splice / merge maybe?
  AtomMapType newAtomMap;

  /* A note on the parallelism:
   * This is safe because we are not actually mutating the unordered_map
   * itself, merely the elements stored in it. The access path towards other
   * elements does not change merely because we mutate the values. For that to
   * happen, we would have to mutate keys, which we aren't doing.
   *
   * Inserting into the new atom map does have to be managed, though.
   */

  auto applyPermutationAndMoveAtomStereopermutator = [&](AtomStereopermutator& permutator) -> void {
    permutator.applyPermutation(permutation);
#pragma omp critical(applyPermutationAtomMapMutation)
    {
      newAtomMap.emplace(permutator.centralIndex(), std::move(permutator));
    }
  };

#pragma omp parallel
  {
#pragma omp single
    {
      for(auto& mapPair : _atomStereopermutators) {
        AtomStereopermutator& permutator = mapPair.second;
#pragma omp task
        {
          applyPermutationAndMoveAtomStereopermutator(permutator);
        }
      }
    }
  }
  std::swap(newAtomMap, _atomStereopermutators);

  /* Create a new bond map (this is not worth parallelizing, applyPermutation
   * on BondStereopermutators is really cheap so the threads just get in each
   * others way mutating newBondMap)
   */
  BondMapType newBondMap;
  for(auto& mapPair : _bondStereopermutators) {
    BondStereopermutator& permutator = mapPair.second;
    permutator.applyPermutation(permutation);
    newBondMap.emplace(permutator.edge(), std::move(permutator));
  }
  std::swap(newBondMap, _bondStereopermutators);
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
  const AtomEnvironmentComponents componentBitmask
) const {
  if(componentBitmask & AtomEnvironmentComponents::Shapes) {
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

      // Ensure the shapes match
      if(stereopermutator.getShape() != otherStereopermutatorOption->getShape()) {
        return false;
      }

      if(componentBitmask & AtomEnvironmentComponents::Stereopermutations) {
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

      if(componentBitmask & AtomEnvironmentComponents::Stereopermutations) {
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

} // namespace Scine
