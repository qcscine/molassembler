/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/StereopermutatorListImpl.h"

#include "temple/Functional.h"
#include "boost/range/adaptor/map.hpp"

namespace Scine {
namespace molassembler {

AtomStereopermutator& StereopermutatorList::Impl::add(
  AtomStereopermutator stereopermutator
) {
  const AtomIndex i = stereopermutator.centralIndex();

  auto emplaceResult = atomStereopermutators.emplace(
    i,
    std::move(stereopermutator)
  );

  if(!emplaceResult.second) {
    throw std::logic_error("Stereopermutator not added. Another is already at its place");
  }

  return emplaceResult.first->second;
}

BondStereopermutator& StereopermutatorList::Impl::add(
  BondStereopermutator stereopermutator
) {
  const BondIndex edge = stereopermutator.edge();

  auto emplaceResult = bondStereopermutators.emplace(
    edge,
    std::move(stereopermutator)
  );

  if(!emplaceResult.second) {
    throw std::logic_error("Stereopermutator not added. Another is already at its place");
  }

  return emplaceResult.first->second;
}


//! Apply an index mapping to the list of stereopermutators
void StereopermutatorList::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
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
      for(auto& mapPair : atomStereopermutators) {
        AtomStereopermutator& permutator = mapPair.second;
#pragma omp task
        {
          applyPermutationAndMoveAtomStereopermutator(permutator);
        }
      }
    }
  }
  std::swap(newAtomMap, atomStereopermutators);

  /* Create a new bond map (this is not worth parallelizing, applyPermutation
   * on BondStereopermutators is really cheap so the threads just get in each
   * others way mutating newBondMap)
   */
  BondMapType newBondMap;
  for(auto& mapPair : bondStereopermutators) {
    BondStereopermutator& permutator = mapPair.second;
    permutator.applyPermutation(permutation);
    newBondMap.emplace(permutator.edge(), std::move(permutator));
  }
  std::swap(newBondMap, bondStereopermutators);
}

void StereopermutatorList::Impl::clear() {
  atomStereopermutators.clear();
  bondStereopermutators.clear();
}

void StereopermutatorList::Impl::clearBonds() {
  bondStereopermutators.clear();
}

void StereopermutatorList::Impl::propagateVertexRemoval(const AtomIndex removedIndex) {
  // Drop any stereopermutators involving this atom from the list
  auto findIter = atomStereopermutators.find(removedIndex);
  if(findIter != atomStereopermutators.end()) {
    atomStereopermutators.erase(findIter);
  }

  /* Go through all state in the StereopermutatorList and decrement any indices
   * larger than the one being removed
   */
  for(auto& stereopermutators : atomStereopermutators | boost::adaptors::map_values) {
    stereopermutators.propagateVertexRemoval(removedIndex);
  }

  /*for(auto& bondMapPair : bondStereopermutators()) {
    bondMapPair.second.propagateVertexRemoval(removedIndex);
  }*/
}

void StereopermutatorList::Impl::remove(const AtomIndex index) {
  auto findIter = atomStereopermutators.find(index);

  if(findIter != atomStereopermutators.end()) {
    atomStereopermutators.erase(findIter);
  } else {
    throw std::logic_error("No such atom stereopermutator found!");
  }
}

void StereopermutatorList::Impl::remove(const BondIndex& edge) {
  auto findIter = bondStereopermutators.find(edge);

  if(findIter != bondStereopermutators.end()) {
    bondStereopermutators.erase(findIter);
  } else {
    throw std::logic_error("No such bond stereopermutator found!");
  }
}

void StereopermutatorList::Impl::try_remove(const AtomIndex index) {
  auto findIter = atomStereopermutators.find(index);

  if(findIter != atomStereopermutators.end()) {
    atomStereopermutators.erase(findIter);
  }
}

void StereopermutatorList::Impl::try_remove(const BondIndex& edge) {
  auto findIter = bondStereopermutators.find(edge);

  if(findIter != bondStereopermutators.end()) {
    bondStereopermutators.erase(findIter);
  }
}

boost::optional<AtomStereopermutator&> StereopermutatorList::Impl::option(const AtomIndex index) {
  auto findIter = atomStereopermutators.find(index);

  if(findIter != atomStereopermutators.end()) {
    return findIter->second;
  }

  return boost::none;
}

boost::optional<BondStereopermutator&> StereopermutatorList::Impl::option(const BondIndex& edge) {
  auto findIter = bondStereopermutators.find(edge);

  if(findIter != bondStereopermutators.end()) {
    return findIter->second;
  }

  return boost::none;
}

/* Information */
bool StereopermutatorList::Impl::empty() const {
  return atomStereopermutators.empty() && bondStereopermutators.empty();
}

unsigned StereopermutatorList::Impl::A() const {
  return atomStereopermutators.size();
}

unsigned StereopermutatorList::Impl::B() const {
  return bondStereopermutators.size();
}

unsigned StereopermutatorList::Impl::size() const {
  return atomStereopermutators.size() + bondStereopermutators.size();
}

boost::optional<const AtomStereopermutator&> StereopermutatorList::Impl::option(const AtomIndex index) const {
  auto findIter = atomStereopermutators.find(index);

  if(findIter != atomStereopermutators.end()) {
    return findIter->second;
  }

  return boost::none;
}

boost::optional<const BondStereopermutator&> StereopermutatorList::Impl::option(const BondIndex& edge) const {
  auto findIter = bondStereopermutators.find(edge);

  if(findIter != bondStereopermutators.end()) {
    return findIter->second;
  }

  return boost::none;
}

bool StereopermutatorList::Impl::hasZeroAssignmentStereopermutators() const {
  return temple::any_of(
    atomStereopermutators | boost::adaptors::map_values,
    [](const auto& stereopermutator) -> bool {
      return stereopermutator.numAssignments() == 0u;
    }
  ) || temple::any_of(
    bondStereopermutators | boost::adaptors::map_values,
    [](const auto& stereopermutator) -> bool {
      return stereopermutator.numAssignments() == 0u;
    }
  );
}

bool StereopermutatorList::Impl::hasUnassignedStereopermutators() const {
  return temple::any_of(
    atomStereopermutators | boost::adaptors::map_values,
    [](const auto& stereopermutator) -> bool {
      return !stereopermutator.assigned();
    }
  ) || temple::any_of(
    bondStereopermutators | boost::adaptors::map_values,
    [](const auto& stereopermutator) -> bool {
      return !stereopermutator.assigned();
    }
  );
}

bool StereopermutatorList::Impl::compare(
  const Impl& other,
  const AtomEnvironmentComponents componentBitmask
) const {
  if(componentBitmask & AtomEnvironmentComponents::Shapes) {
    // Check sizes
    if(
      atomStereopermutators.size() != other.atomStereopermutators.size()
      || bondStereopermutators.size() != other.bondStereopermutators.size()
    ) {
      return false;
    }

    // Check all atom stereopermutators
    for(const auto& stereopermutator : atomStereopermutators | boost::adaptors::map_values) {
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
    for(const auto& stereopermutator : bondStereopermutators | boost::adaptors::map_values) {
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
bool StereopermutatorList::Impl::operator == (const Impl& other) const {
  return (
    atomStereopermutators == other.atomStereopermutators
    && bondStereopermutators == other.bondStereopermutators
  );
}

} // namespace molassembler
} // namespace Scine
