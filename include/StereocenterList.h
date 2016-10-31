#ifndef INCLUDE_GRAPH_FEATURE_LIST_H
#define INCLUDE_GRAPH_FEATURE_LIST_H

#include <map>

#include "Stereocenter.h"

namespace MoleculeManip {

class StereocenterList {
private:
/* Private members */
  std::vector<
    std::shared_ptr<
      Stereocenters::Stereocenter
    >
  > _features;

public:
  /* Modification */
  /*!
   * Removes all stored Stereocenters
   */
  void invalidate() {
    _features.clear();
  }

  /*!
   * Removes Stereocenters centered on a specific index
   */
  void selectivelyInvalidate(const AtomIndexType& a) {
    _features.erase(
      std::remove_if(
        _features.begin(),
        _features.end(),
        [&a](const std::shared_ptr<Stereocenters::Stereocenter>& featurePtr) {
          return (
            featurePtr -> involvedAtoms()
          ).count(a) == 1;
        }
      ),
      _features.end()
    );
  }

  /*!
   * Adds a shared_ptr instance to the list.
   */
  void add(const std::shared_ptr<Stereocenters::Stereocenter>& featurePtr) {
    _features.push_back(featurePtr);
  }

  /* Information */
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

  /* Iterators */
  // Begin and end const iterators for easy traversal
  std::vector<
    std::shared_ptr<
      Stereocenters::Stereocenter
    >
  >::const_iterator begin() const {
    return _features.cbegin();
  }
  std::vector<
    std::shared_ptr<
      Stereocenters::Stereocenter
    >
  >::const_iterator end() const {
    return _features.cend();
  }

};

}

#endif
