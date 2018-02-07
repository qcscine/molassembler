#include "GenerateUniques.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <sstream>

namespace UniqueAssignments {

bool predicateHasTransArrangedPairs(
  const Assignment& assignment,
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


std::vector<Assignment> uniqueAssignments(
  const Assignment& initial,
  const Symmetry::Name& symmetryName,
  const bool& removeTransSpanningGroups
) {
  /* NOTE: This algorithm may seem wasteful in terms of memory (after all, one
   * could, insted of keeping a full set of all rotations of all unique 
   * assignments, just do pair-wise comparisons between a new assignment and 
   * all existing unique ones. However, one would be doing a lot of repeated 
   * work, since the pair-wise comparison (see isRotationallySuperimposable) 
   * just generates rotations of one and compares those with the other. It is 
   * chosen here to prefer speed over memory requirements. After all, the 
   * number of Assignment objects that will be generated and stored is unlikely
   * to pass 1000.
   */

  // make a copy of initial so we can modify it by permutation
  Assignment assignment = initial;

  // ensure we start with the lowest permutation
  assignment.lowestPermutation();

  /* in case we want to skip trans pairs, the initial assignment must also not
   * have any trans spanning pairs
   */
  if(removeTransSpanningGroups) {
    while(predicateHasTransArrangedPairs(assignment, symmetryName)) {
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
  std::vector<Assignment> uniqueAssignments {*rotationsSet.begin()};

  // go through all possible permutations of columns
  while(assignment.nextPermutation()) {
    if( // skip permutations with trans pairs if desired
      removeTransSpanningGroups 
      && predicateHasTransArrangedPairs(assignment, symmetryName)
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
      uniqueAssignments.push_back(*assignmentRotations.begin());

      // and add its rotations to the set
      rotationsSet.insert(
        assignmentRotations.begin(),
        assignmentRotations.end()
      );

    } 
  }
    
  return uniqueAssignments;
}

UniqueAssignmentsReturnType uniqueAssignmentsWithCounts(
  const Assignment& initial,
  const Symmetry::Name& symmetryName,
  const bool& removeTransSpanningGroups
) {
  /* NOTE: This algorithm may seem wasteful in terms of memory (after all, one
   * could, insted of keeping a full set of all rotations of all unique 
   * assignments, just do pair-wise comparisons between a new assignment and 
   * all existing unique ones. However, one would be doing a lot of repeated 
   * work, since the pair-wise comparison (see isRotationallySuperimposable) 
   * just generates rotations of one and compares those with the other. It is 
   * here chosen to prefer speed over memory requirements. After all, the 
   * number of Assignment objects that will be generated and stored is unlikely
   * to pass 1000.
   */

  struct RotationsAndOccurrencesCount {
    std::set<Assignment> rotations;
    unsigned occurrencesCount = 1;

    // Help constructor
    RotationsAndOccurrencesCount() = default;
    explicit RotationsAndOccurrencesCount(std::set<Assignment>&& rotations) : rotations(rotations)
    {}
  };

  // make a copy of initial so we can modify it by permutation
  Assignment assignment = initial;

  // ensure we start with the lowest permutation
  assignment.lowestPermutation();

  /* in case we want to skip trans pairs, the initial assignment must also not
   * have any trans spanning pairs
   */
  if(removeTransSpanningGroups) {
    while(predicateHasTransArrangedPairs(assignment, symmetryName)) {
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
  Assignment lowestRotation = *initialRotations.begin();

  std::map<
    Assignment,
    RotationsAndOccurrencesCount
  > trackingData {
    {
      lowestRotation,
      RotationsAndOccurrencesCount {std::move(initialRotations)}
    }
  };

  // go through all possible permutations of columns
  while(assignment.nextPermutation()) {
    if( // skip permutations with trans pairs if desired
      removeTransSpanningGroups 
      && predicateHasTransArrangedPairs(assignment, symmetryName)
    ) {
      continue;
    }

    // is the current assignment not contained within the set of rotations?
    bool isContained = false;
    for(auto& iterPair : trackingData) {
      if(iterPair.second.rotations.count(assignment) != 0) {
        isContained = true;
        iterPair.second.occurrencesCount += 1;
      }
    }

    if(!isContained) {
      auto rotations = assignment.generateAllRotations(symmetryName);
      auto lowestRotation = *rotations.begin();
      trackingData[lowestRotation] = RotationsAndOccurrencesCount {
        std::move(rotations)
      };
    }

  }

  // Condense the data down to vectors
  UniqueAssignmentsReturnType returnData;

  // Get the data
  for(auto& iterPair: trackingData) {
    returnData.assignments.push_back(
      iterPair.first
    );
    returnData.weights.push_back(
      iterPair.second.occurrencesCount
    );
  }

  return returnData;
}

} // namespace UniqueAssignments
