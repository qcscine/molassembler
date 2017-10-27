#ifndef INCLUDE_GRAPH_FEATURE_LIST_H
#define INCLUDE_GRAPH_FEATURE_LIST_H

#include <map>

#include "CNStereocenter.h"
#include "EZStereocenter.h"

/*! @file
 *
 * Contains the declaration for a class that stores a list of all stereocenters
 * in a molecule.
 */

namespace MoleculeManip {

class StereocenterMap {
public:
  using PtrType = std::shared_ptr<Stereocenters::Stereocenter>;
  using MapType = std::map<
    AtomIndexType,
    PtrType
  >;

private:
  MapType _stereocenters;

public:
/* Modifiers */
  StereocenterMap() = default;
  StereocenterMap(std::initializer_list<PtrType> initList) {
    for(const auto& ptr : initList) {
      add(ptr);
    }
  }

  //! Add a stereocenter
  void add(const PtrType& stereocenterPtr) {
    // Ensure the involved atoms aren't mapped yet
    for(const auto& atomIndex : stereocenterPtr -> involvedAtoms()) {
      if(_stereocenters.count(atomIndex) == 1) {
        throw std::logic_error(
          "Trying to add a stereocenter to an atom that is already mapped!"
        );
      }
    }

    for(const auto& atomIndex : stereocenterPtr -> involvedAtoms()) {
      _stereocenters[atomIndex] = stereocenterPtr;
    }
  }

  //! Removes all stereocenters
  void clear() noexcept {
    _stereocenters.clear();
  }

  /*! 
   * Removes a pointer if it is the SAME as the one passed (using address),
   * does not consider equality of the types involved
   */
  void remove(const PtrType& stereocenterPtr) {
    for(const auto& index : stereocenterPtr -> involvedAtoms()) {
      assert(
        _stereocenters.count(index) == 1
        && _stereocenters.at(index) == stereocenterPtr
      );

      _stereocenters.erase(index);
    }
  }

/* Information */
  bool empty() const noexcept {
    return _stereocenters.empty();
  }

  //! Returns an optional with the stereocenter pointer involving the atom index
  boost::optional<PtrType> get(const AtomIndexType& index) const {
    if(_stereocenters.count(index) == 1) {
      return _stereocenters.at(index);
    }

    return boost::none;
  }

  unsigned size() const noexcept {
    return _stereocenters.size();
  }

  MapType::iterator begin() {
    return _stereocenters.begin();
  }

  MapType::iterator end() {
    return _stereocenters.end();
  }

  MapType::const_iterator cbegin() const {
    return _stereocenters.cbegin();
  }

  MapType::const_iterator cend() const {
    return _stereocenters.cend();
  }

  bool operator == (const StereocenterMap& other) const {
    return std::equal(
      _stereocenters,
      other._stereocenters,
      [&](const auto& iterPairA, const auto& iterPairB) -> bool {
        return (
          iterPairA.first == iterPairB.first
          && Stereocenters::compareStereocenterEqual(
            iterPairA.second,
            iterPairB.second
          )
        );
      }
    );
  }
};

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
            return Stereocenters::compareStereocenterEqual(stereocenter, otherStereocenter);
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
