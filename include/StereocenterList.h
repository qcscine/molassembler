#ifndef INCLUDE_GRAPH_FEATURE_LIST_H
#define INCLUDE_GRAPH_FEATURE_LIST_H

#include <map>

#include "Stereocenter.h"
#include "CNStereocenter.h"
#include "EZStereocenter.h"

namespace MoleculeManip {

class StereocenterList {
public:
  using PtrType = std::shared_ptr<Stereocenters::Stereocenter>;
  /* Important detail here:
   * Remember that when using shared_ptr in an ordered container, all comparison
   * operators simply compare the address where elements are stored. When 
   * checking equality, this merely amounts to checking if they are the SAME, 
   * NOT whether they are EQUAL.
   */
  using SetType = std::set<PtrType>;

  using ListType = std::vector<PtrType>;

  using const_iterator = SetType::const_iterator;
  using iterator = SetType::iterator;

private:
/* Private members */
  SetType _features;

public:
/* Constructors */
  StereocenterList() = default;
  StereocenterList(std::initializer_list<PtrType> initList) {
    for(const auto& ptr : initList) {
      add(ptr);
    }
  }

  /* Modification */
  /*!
   * Removes all stored Stereocenters
   */
  void invalidate() {
    _features.clear();
  }

  /*!
   * Adds a shared_ptr instance to the list.
   */
  void add(const std::shared_ptr<Stereocenters::Stereocenter>& featurePtr) {
    _features.insert(featurePtr);
  }

  /* Information */
  boost::optional<
    std::shared_ptr<Stereocenters::Stereocenter>
  > involving(const AtomIndexType& index) {
    for(const auto& stereocenterPtr : _features) {
      if(stereocenterPtr -> involvedAtoms().count(index)) {
        return stereocenterPtr;
      }
    }

    return boost::none;
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
   * Returns a map of atom indices to vectors of Stereocenters. 
   * The returned map does NOT contain a key for every atom in the Molecule.
   * When using the returned map, use count() before accessing with at()
   */
  std::map<
    AtomIndexType,
    std::vector<
      std::shared_ptr<
        Stereocenters::Stereocenter
      >
    >
  > getAtomIndexMapping() const {
    std::map<
      AtomIndexType,
      std::vector<
        std::shared_ptr<
          Stereocenters::Stereocenter
        >
      >
    > mapping;

    for(const auto& featurePtr : _features) {
      auto atomSet = featurePtr -> involvedAtoms();

      for(const auto& atom : atomSet) {
        if(mapping.count(atom) == 0) {
          mapping[atom] = {featurePtr};
        } else {
          mapping[atom].push_back(featurePtr);
        }
      }

    }

    return mapping;
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

  bool operator == (const StereocenterList& other) {
    using namespace Stereocenters;

    /* TODO
     * - This is inefficient because items that have been matched to from a in
     *   b are still considered when the next item in a is sought:
     *
     *   a  b
     *   1  1
     *   2  2
     *
     *   1. look for a's 1 in b -> look at 1, found
     *   2. look for a's 2 in b -> *look at 1*, look at 2, found
     */
    if(_features.size() != other._features.size()) return false;

    bool all_found = true;
    for(const auto& stereocenter: _features) {
      if(
        std::find_if(
          other._features.begin(),
          other._features.end(),
          [&](const std::shared_ptr<Stereocenter>& otherStereocenter) -> bool {
            return Stereocenters::strictComparePtr(stereocenter, otherStereocenter);
          }
        ) == other._features.end()
      ) {
        all_found = false;
        break;
      }
    }

    return all_found;

  }
};

}

#endif
