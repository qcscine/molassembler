#ifndef INCLUDE_GRAPH_FEATURE_LIST_H
#define INCLUDE_GRAPH_FEATURE_LIST_H

#include <map>

#include "Stereocenter.h"
#include "CNStereocenter.h"
#include "EZStereocenter.h"

namespace MoleculeManip {

class StereocenterList {
public:
  using ListType = std::vector<
    std::shared_ptr<
      Stereocenters::Stereocenter
    >
  >;

private:
/* Private members */
  ListType _features;

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

  unsigned size() const {
    return _features.size();
  }

  /* Iterators */
  // Begin and end iterators for easy traversal
  ListType::iterator begin() {
    return _features.begin();
  }

  ListType::iterator end() {
    return _features.end();
  }

  ListType::const_iterator begin() const {
    return _features.cbegin();
  }

  ListType::const_iterator end() const {
    return _features.cend();
  }

  bool operator == (const StereocenterList& other) {
    using namespace Stereocenters;

    auto stereocentersEqual = [](
      const std::shared_ptr<Stereocenter>& a,
      const std::shared_ptr<Stereocenter>& b
    ) -> bool {
      if(a -> type() == b -> type()) {
        if(a -> type() == Type::CNStereocenter) {
          auto aDerived = std::dynamic_pointer_cast<CNStereocenter>(a);
          auto bDerived = std::dynamic_pointer_cast<CNStereocenter>(b);

          return *aDerived == *bDerived;
        } else { // EZStereocenter
          auto aDerived = std::dynamic_pointer_cast<EZStereocenter>(a);
          auto bDerived = std::dynamic_pointer_cast<EZStereocenter>(b);

          return *aDerived == *bDerived;
        }
      } else {
        return false;
      }
    };

    if(_features.size() != other._features.size()) return false;

    bool all_found = true;
    for(const auto& stereocenter: _features) {
      if(
        std::find_if(
          other._features.begin(),
          other._features.end(),
          [&](const std::shared_ptr<Stereocenter>& otherStereocenter) -> bool {
            return stereocentersEqual(stereocenter, otherStereocenter);
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
