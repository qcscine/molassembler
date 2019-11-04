/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "stereopermutation/GenerateUniques.h"

#include "boost/functional/hash.hpp"

#include "shapes/Data.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <sstream>
#include <unordered_map>

namespace Scine {

namespace stereopermutation {

bool hasTransArrangedPairs(
  const Stereopermutation& stereopermutation,
  const Symmetry::Shape shape
) {
  // for every pair in links
  for(const auto& indexPair : stereopermutation.links) {
    if(
      Symmetry::angleFunction(shape)(
        indexPair.first,
        indexPair.second
      ) == M_PI
    ) {
      return true;
    }
  }

  return false;
}


std::vector<Stereopermutation> uniques(
  const Stereopermutation& initial,
  const Symmetry::Shape shape,
  const bool removeTransSpanningGroups
) {
  /* NOTE: This algorithm may seem wasteful in terms of memory (after all, one
   * could, insted of keeping a full set of all rotations of all unique
   * stereopermutations, just do pair-wise comparisons between a new stereopermutation and
   * all existing unique ones. However, one would be doing a lot of repeated
   * work, since the pair-wise comparison (see isRotationallySuperimposable)
   * just generates rotations of one and compares those with the other. It is
   * chosen here to prefer speed over memory requirements. After all, the
   * number of Stereopermutation objects that will be generated and stored is unlikely
   * to pass 1000.
   */

  // make a copy of initial so we can modify it by permutation
  Stereopermutation stereopermutation = initial;

  // ensure we start with the lowest permutation
  stereopermutation.lowestPermutation();

  /* in case we want to skip trans pairs, the initial stereopermutation must also not
   * have any trans spanning pairs
   */
  if(removeTransSpanningGroups) {
    while(hasTransArrangedPairs(stereopermutation, shape)) {
      bool hasAnotherPermutation = stereopermutation.nextPermutation();
      if(!hasAnotherPermutation) {
        /* This can happen, e.g. in square-planar AAAB with
         * links: {0, 3}, {1, 3}, {2, 3}, every possible permutation contains
         * trans-arranged pairs. Then we return an empty vector.
         */
        return {};
      }
    }
  }

  // Generate all rotations of the initial stereopermutation
  auto rotationsSet = stereopermutation.generateAllRotations(shape);

  // The lowest rotation of the passed stereopermutation is the first unique stereopermutation
  std::vector<Stereopermutation> uniqueStereopermutations {*rotationsSet.begin()};

  // go through all possible permutations of columns
  while(stereopermutation.nextPermutation()) {
    if( // skip permutations with trans pairs if desired
      removeTransSpanningGroups
      && hasTransArrangedPairs(stereopermutation, shape)
    ) {
      continue;
    }

    // is the current stereopermutation not contained within the set of rotations?
    if(
      rotationsSet.count(stereopermutation) == 0
    ) {
      // if so, it is a unique stereopermutation, generate all rotations
      auto stereopermutationRotations = stereopermutation.generateAllRotations(shape);

      // add the smallest stereopermutation from the generated set to the list of uniques
      uniqueStereopermutations.push_back(*stereopermutationRotations.begin());

      // and add its rotations to the set
      rotationsSet.insert(
        stereopermutationRotations.begin(),
        stereopermutationRotations.end()
      );

    }
  }

  return uniqueStereopermutations;
}

StereopermutationsWithWeights uniquesWithWeights(
  const Stereopermutation& initial,
  const Symmetry::Shape shape,
  const bool removeTransSpanningGroups
) {
  /* NOTE: This algorithm may seem wasteful in terms of memory (after all, one
   * could, insted of keeping a full set of all rotations of all unique
   * stereopermutations, just do pair-wise comparisons between a new stereopermutation and
   * all existing unique ones. However, one would be doing a lot of repeated
   * work, since the pair-wise comparison (see isRotationallySuperimposable)
   * just generates rotations of one and compares those with the other. It is
   * here chosen to prefer speed over memory requirements. After all, the
   * number of Stereopermutation objects that will be generated and stored is unlikely
   * to pass 1000.
   */

  struct RotationsAndOccurrencesCount {
    std::set<Stereopermutation> rotations;
    unsigned occurrencesCount = 1;

    // Help constructor
    RotationsAndOccurrencesCount() = default;
    explicit RotationsAndOccurrencesCount(std::set<Stereopermutation> passRotations)
      : rotations(std::move(passRotations)) {}
  };

  // make a copy of initial so we can modify it by permutation
  Stereopermutation stereopermutation = initial;

  // ensure we start with the lowest permutation
  stereopermutation.lowestPermutation();

  /* in case we want to skip trans pairs, the initial stereopermutation must also not
   * have any trans spanning pairs
   */
  if(removeTransSpanningGroups) {
    while(hasTransArrangedPairs(stereopermutation, shape)) {
      bool hasAnotherPermutation = stereopermutation.nextPermutation();
      if(!hasAnotherPermutation) {
        /* This can happen, e.g. in square-planar AAAB with
         * links: {0, 3}, {1, 3}, {2, 3}, every possible permutation contains
         * trans-arranged pairs. Then we return an empty vector.
         */
        return {};
      }
    }
  }

  auto initialRotations = stereopermutation.generateAllRotations(shape);
  const Stereopermutation& lowestRotation = *initialRotations.begin();

  StereopermutationsWithWeights data;

  data.weights.reserve(40);
  data.stereopermutations.reserve(40);

  data.stereopermutations.push_back(lowestRotation);
  data.weights.push_back(1);

  std::unordered_map<Stereopermutation, unsigned, boost::hash<Stereopermutation>> rotationCounterMap;

  for(const auto& rotation: initialRotations) {
    rotationCounterMap.emplace(
      rotation,
      0
    );
  }

  while(stereopermutation.nextPermutation()) {
    if( // skip permutations with trans pairs if desired
      removeTransSpanningGroups
      && hasTransArrangedPairs(stereopermutation, shape)
    ) {
      continue;
    }

    auto findIter = rotationCounterMap.find(stereopermutation);
    if(findIter == rotationCounterMap.end()) {
      auto rotations = stereopermutation.generateAllRotations(shape);
      data.stereopermutations.push_back(*rotations.begin());
      data.weights.push_back(1);

      for(const auto& rotation : rotations) {
        rotationCounterMap.emplace(
          rotation,
          data.weights.size() - 1
        );
      }
    } else {
      ++data.weights.at(findIter->second);
    }
  }

  data.weights.shrink_to_fit();
  data.stereopermutations.shrink_to_fit();

  return data;
}

} // namespace stereopermutation

} // namespace Scine
