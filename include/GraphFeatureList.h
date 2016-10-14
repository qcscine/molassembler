#ifndef INCLUDE_GRAPH_FEATURE_LIST_H
#define INCLUDE_GRAPH_FEATURE_LIST_H

#include <map>

#include "GraphFeature.h"

namespace MoleculeManip {

class GraphFeatureList {
private:
  std::vector<
    std::shared_ptr<
      GraphFeatures::GraphFeature
    >
  > _features;

public:
  /* Modification */
  void invalidate() {
    _features.clear();
  }

  void selectivelyInvalidate(const AtomIndexType& a) {
    _features.erase(
      std::remove_if(
        _features.begin(),
        _features.end(),
        [&a](const std::shared_ptr<GraphFeatures::GraphFeature>& featurePtr) {
          return (
            featurePtr -> involvedAtoms()
          ).count(a) == 1;
        }
      ),
      _features.end()
    );
  }

  void add(const std::shared_ptr<GraphFeatures::GraphFeature>& featurePtr) {
    _features.push_back(featurePtr);
  }

  /* Information */
  /*!
   * The returned map does NOT contain a key for every atom in the Molecule.
   * When using the returned map, use count() before accessing with at()
   */
  std::map<
    AtomIndexType,
    std::vector<
      std::shared_ptr<
        GraphFeatures::GraphFeature
      >
    >
  > getAtomIndexMapping() const {
    std::map<
      AtomIndexType,
      std::vector<
        std::shared_ptr<
          GraphFeatures::GraphFeature
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
  std::vector<
    std::shared_ptr<
      GraphFeatures::GraphFeature
    >
  >::const_iterator begin() const {
    return _features.cbegin();
  }
  std::vector<
    std::shared_ptr<
      GraphFeatures::GraphFeature
    >
  >::const_iterator end() const {
    return _features.cend();
  }

};

}

#endif
