#ifndef LIB_UNIQUE_ASSIGNMENTS_GENERATE_UNIQUES_H
#define LIB_UNIQUE_ASSIGNMENTS_GENERATE_UNIQUES_H

#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
#include <sstream>

#include "Assignment.h"

/* TODO
 */

namespace UniqueAssignments {

template<class Symmetry>
bool predicateHasTransArrangedPairs(
  const Assignment<Symmetry>& assignment
) {
  // for every pair in links
  for(const auto& indexPair : assignment.links) {
    if(
      Symmetry::angle(
        indexPair.first,
        indexPair.second
      ) == 180.0
    ) return true;
  }
  return false;
}

/* NOTE: Gives NO guarantees as to satisfiability (if assignments can be
 *  fulfilled with real ligands) 
 * E.g. M (A-A)_3 generates a trans-trans-trans assignment, which is extremely 
 *  hard to find actual ligands for that work.
 * The satisfiability of assignments must be checked before trying to embed 
 *  structures with completely nonsensical constraints. Perhaps restrict A-A 
 *  ligands with bridge length 4 (chelating atoms included), maybe even up to 6
 *  to cis arrangements. Xantphos (with bridge length 7) is the smallest 
 *  trans-spanning ligand mentioned in Wikipedia.
 */
template<class Symmetry>
std::vector<
  Assignment<Symmetry>
> uniqueAssignments(
  const Assignment<Symmetry>& initial,
  const bool& removeTransSpanningGroups = true
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

  /* TODO: 
   * - If the rotationsSet size reaches 720, isn't it impossible for there to 
   *   be more uniques? Can we skip the remaining permutations?
   */

  // make a copy of initial so we can modify it by permutation
  Assignment<Symmetry> assignment = initial;

  // ensure we start with the lowest permutation
  assignment.lowestPermutation();

  /* in case we want to skip trans pairs, the initial assignment must also not
   * have any trans spanning pairs
   */
  if(removeTransSpanningGroups) {
    while(predicateHasTransArrangedPairs<Symmetry>(assignment)) {
      assignment.nextPermutation();
    }
  }

  // The provided assignment is the first unique assignment
  std::vector<
    Assignment<Symmetry>
  > uniqueAssignments {assignment};

  // generate the initial assignment's set of rotations
  auto rotationsSet = assignment.generateAllRotations();

  // go through all possible permutations of columns
  while(assignment.nextPermutation()) {
    if(
      removeTransSpanningGroups 
      && predicateHasTransArrangedPairs<Symmetry>(assignment)
    ) continue; // skip permutations with trans pairs if desired

    // is the current assignment not contained within the set of rotations?
    if(
      rotationsSet.count(assignment) == 0
    ) {
      // if so, it is a unique assignment, so add it to the list
      uniqueAssignments.push_back(assignment);

      // and add its rotations to the set
      /* C++17
      rotationsSet.merge(
        assignment.generateAllRotations()
      ); */
      auto assignmentRotations = assignment.generateAllRotations();
      rotationsSet.insert(
        assignmentRotations.begin(),
        assignmentRotations.end()
      );

    } 
  }
    
  return uniqueAssignments;
}

} // eo namespace

#endif
