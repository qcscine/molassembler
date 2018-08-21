#include "stereopermutation/GenerateUniques.h"

#include "boost/functional/hash.hpp"

#include "chemical_symmetries/Symmetries.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <sstream>
#include <unordered_map>

namespace stereopermutation {

bool hasTransArrangedPairs(
  const Stereopermutation& assignment,
  const Symmetry::Name& symmetryName
) {
  // for every pair in links
  for(const auto& indexPair : assignment.links) {
    if(
      Symmetry::angleFunction(symmetryName)(
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
  const Symmetry::Name& symmetryName,
  const bool removeTransSpanningGroups
) {
  /* NOTE: This algorithm may seem wasteful in terms of memory (after all, one
   * could, insted of keeping a full set of all rotations of all unique
   * assignments, just do pair-wise comparisons between a new assignment and
   * all existing unique ones. However, one would be doing a lot of repeated
   * work, since the pair-wise comparison (see isRotationallySuperimposable)
   * just generates rotations of one and compares those with the other. It is
   * chosen here to prefer speed over memory requirements. After all, the
   * number of Stereopermutation objects that will be generated and stored is unlikely
   * to pass 1000.
   */

  // make a copy of initial so we can modify it by permutation
  Stereopermutation assignment = initial;

  // ensure we start with the lowest permutation
  assignment.lowestPermutation();

  /* in case we want to skip trans pairs, the initial assignment must also not
   * have any trans spanning pairs
   */
  if(removeTransSpanningGroups) {
    while(hasTransArrangedPairs(assignment, symmetryName)) {
      bool hasAnotherPermutation = assignment.nextPermutation();
      if(!hasAnotherPermutation) {
        /* This can happen, e.g. in square-planar AAAB with
         * links: {0, 3}, {1, 3}, {2, 3}, every possible permutation contains
         * trans-arranged pairs. Then we return an empty vector.
         */
        return {};
      }
    }
  }

  // Generate all rotations of the initial assignment
  auto rotationsSet = assignment.generateAllRotations(symmetryName);

  // The lowest rotation of the passed assignment is the first unique assignment
  std::vector<Stereopermutation> uniqueStereopermutations {*rotationsSet.begin()};

  // go through all possible permutations of columns
  while(assignment.nextPermutation()) {
    if( // skip permutations with trans pairs if desired
      removeTransSpanningGroups
      && hasTransArrangedPairs(assignment, symmetryName)
    ) {
      continue;
    }

    // is the current assignment not contained within the set of rotations?
    if(
      rotationsSet.count(assignment) == 0
    ) {
      // if so, it is a unique assignment, generate all rotations
      auto assignmentRotations = assignment.generateAllRotations(symmetryName);

      // add the smallest assignment from the generated set to the list of uniques
      uniqueStereopermutations.push_back(*assignmentRotations.begin());

      // and add its rotations to the set
      rotationsSet.insert(
        assignmentRotations.begin(),
        assignmentRotations.end()
      );

    }
  }

  return uniqueStereopermutations;
}

StereopermutationsWithWeights uniquesWithWeights(
  const Stereopermutation& initial,
  const Symmetry::Name& symmetryName,
  const bool removeTransSpanningGroups
) {
  /* NOTE: This algorithm may seem wasteful in terms of memory (after all, one
   * could, insted of keeping a full set of all rotations of all unique
   * assignments, just do pair-wise comparisons between a new assignment and
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
    explicit RotationsAndOccurrencesCount(std::set<Stereopermutation>&& rotations) : rotations(rotations) {}
  };

  // make a copy of initial so we can modify it by permutation
  Stereopermutation assignment = initial;

  // ensure we start with the lowest permutation
  assignment.lowestPermutation();

  /* in case we want to skip trans pairs, the initial assignment must also not
   * have any trans spanning pairs
   */
  if(removeTransSpanningGroups) {
    while(hasTransArrangedPairs(assignment, symmetryName)) {
      bool hasAnotherPermutation = assignment.nextPermutation();
      if(!hasAnotherPermutation) {
        /* This can happen, e.g. in square-planar AAAB with
         * links: {0, 3}, {1, 3}, {2, 3}, every possible permutation contains
         * trans-arranged pairs. Then we return an empty vector.
         */
        return {};
      }
    }
  }

  auto initialRotations = assignment.generateAllRotations(symmetryName);
  const Stereopermutation& lowestRotation = *initialRotations.begin();

  StereopermutationsWithWeights data;

  data.weights.reserve(40);
  data.assignments.reserve(40);

  data.assignments.push_back(lowestRotation);
  data.weights.push_back(1);

  std::unordered_map<Stereopermutation, unsigned, boost::hash<Stereopermutation>> rotationCounterMap;

  for(const auto& rotation: initialRotations) {
    rotationCounterMap.emplace(
      rotation,
      0
    );
  }

  while(assignment.nextPermutation()) {
    if( // skip permutations with trans pairs if desired
      removeTransSpanningGroups
      && hasTransArrangedPairs(assignment, symmetryName)
    ) {
      continue;
    }

    auto findIter = rotationCounterMap.find(assignment);
    if(findIter == rotationCounterMap.end()) {
      auto rotations = assignment.generateAllRotations(symmetryName);
      data.assignments.push_back(*rotations.begin());
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
  data.assignments.shrink_to_fit();

  return data;
}

} // namespace stereopermutation
