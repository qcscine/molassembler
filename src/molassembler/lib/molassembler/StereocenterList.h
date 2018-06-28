#ifndef INCLUDE_GRAPH_FEATURE_LIST_H
#define INCLUDE_GRAPH_FEATURE_LIST_H

#include <unordered_map>

#include "AtomStereocenter.h"
#include "BondStereocenter.h"

/*! @file
 *
 * Contains the declaration for a class that stores a list of all stereocenters
 * in a molecule.
 */

namespace molassembler {

class StereocenterList {
public:
  using PtrType = std::shared_ptr<Stereocenters::Stereocenter>;
  /* Important detail here:
   * Remember that when using shared_ptr in an ordered container, all comparison
   * operators simply compare the address where elements are stored. When
   * checking equality, this merely amounts to checking if they are the SAME,
   * NOT whether they are EQUAL.
   */

  struct PtrCompare {
    bool operator () (const PtrType& a, const PtrType& b) const {
      return Stereocenters::compareStereocenterLessThan(a, b);
    }
  };

  using SetType = std::set<PtrType, PtrCompare>;
  using ListType = std::vector<PtrType>;
  using MapType = std::unordered_map<AtomIndexType, PtrType>;

  using const_iterator = SetType::const_iterator;
  using iterator = SetType::iterator;

private:
/* Private members */
  SetType _features;
  MapType _indexMap;

public:
/* Constructors */
  StereocenterList() = default;
  StereocenterList(std::initializer_list<PtrType> initList) {
    for(const auto& ptr : initList) {
      add(ptr);
    }
  }

/* Modification */
  //! Removes all stored Stereocenters
  void clear() {
    _features.clear();
    _indexMap.clear();
  }

  //! Adds a shared_ptr instance to the list.
  void add(const std::shared_ptr<Stereocenters::Stereocenter>& featurePtr) {
    _features.insert(featurePtr);

    // Add to map as well
    for(const auto& involvedAtom : featurePtr -> involvedAtoms()) {
      _indexMap.emplace(
        involvedAtom,
        featurePtr
      );
    }
  }

  PtrType at(const AtomIndexType index) {
    return _indexMap.at(index);
  }

  /* Removing a vertex invalidates some vertex descriptors, which are used
   * liberally in the stereocenter classes. This function ensures that
   * vertex descriptors are valid throughout.
   */
  void propagateVertexRemoval(const AtomIndexType removedIndex) {
    // Drop any stereocenters involving this atom from the list
    if(involving(removedIndex)) {
      remove(removedIndex);
    }

    /* Go through all state in the StereocenterList and decrement any indices
     * larger than the one being removed
     */
    for(const auto& stereocenterPtr : _features) {
      stereocenterPtr -> propagateVertexRemoval(removedIndex);
    }

    // Rebuild the index map
    _indexMap.clear();
    for(const auto& stereocenterPtr : _features) {
      for(const auto& involvedAtom : stereocenterPtr -> involvedAtoms()) {
        _indexMap.emplace(
          involvedAtom,
          stereocenterPtr
        );
      }
    }
  }

  void remove(const AtomIndexType index) {
    if(_indexMap.count(index) == 0) {
      throw std::logic_error("StereocenterList::remove: No mapping for that index!");
    }

    _features.erase(_indexMap.at(index));
    _indexMap.erase(index);
  }

/* Information */
  const PtrType at(const AtomIndexType index) const {
    return _indexMap.at(index);
  }

  bool involving(const AtomIndexType index) const {
    if(_indexMap.count(index) > 0) {
      return true;
    }

    return false;
  }

  bool isStereogenic(const AtomIndexType index) const {
    auto findIter = _indexMap.find(index);
    if(findIter == _indexMap.end()) {
      return false;
    }

    return findIter->second->numStereopermutations() > 1;
  }
  /* As long as the constraint exists that Stereocenters' involvedAtoms do not
   * overlap, this variant below is unneeded.
   */
  /*ListType involving(const AtomIndexType& index) {
    ListType matches;

    for(const auto& stereocenterPtr : _features) {
      if(stereocenterPtr -> involvedAtoms().count(index)) {
        matches.emplace_back(stereocenterPtr);
      }
    }

    return matches;
  }*/

  /*!
   * Returns a map of atom indices to Stereocenter pointers.
   * The returned map does NOT contain a key for every atom in the Molecule.
   * When using the returned map, use count() before accessing with at()
   */
  const MapType& getAtomIndexMapping() const {
    return _indexMap;
  }

  bool empty() const {
    return _features.empty();
  }

  unsigned size() const {
    return _features.size();
  }

  /* Iterators */
  // Begin and end iterators for easy traversal
  SetType::iterator begin() {
    return _features.begin();
  }

  SetType::iterator end() {
    return _features.end();
  }

  SetType::const_iterator begin() const {
    return _features.cbegin();
  }

  SetType::const_iterator end() const {
    return _features.cend();
  }

  bool operator == (const StereocenterList& other) const {
    return std::equal(
      _features.begin(),
      _features.end(),
      other._features.begin(),
      other._features.end(),
      [](const auto& ptrA, const auto& ptrB) -> bool {
        return Stereocenters::compareStereocenterEqual(ptrA, ptrB);
      }
    );
  }
};

}

#endif
