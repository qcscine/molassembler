#ifndef LIB_UNIQUE_ASSIGNMENTS_GENERATE_UNIQUES_H
#define LIB_UNIQUE_ASSIGNMENTS_GENERATE_UNIQUES_H

#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
#include <iostream>
#include <sstream>

#include "Assignment.h"

/* TODO
 */

namespace UniqueAssignments {

template<
  template<typename T = AssignmentColumn>
  class Symmetry
>
bool predicateHasTransArrangedPairs(
  const Assignment<Symmetry>& assignment
) {
  // for every group in the positionOccupations
  for(
    unsigned i = 0;
    i < assignment.positionOccupations[0].groups.size();
    i++
  ) {
    std::vector<unsigned> pair;
    for(unsigned j = 0; j < Symmetry<>::size; j++) {
      // i is the row index, j the column index
      if(assignment.positionOccupations.at(j).groups.at(i)) pair.push_back(j);
    }
    if(
      Symmetry<>::angle(
        pair[0],
        pair[1]
      ) == 180.0
    ) {
      return true;
    }
  }

  return false;
}

/* Gives NO guarantees as to satisfiability (if assignments can be fulfilled 
 * with real ligands) 
 * E.g. M (A-A)_3 generates a trans-trans-trans assignment, which is extremely 
 *  hard to find actual ligands for that work.
 * The satisfiability of assignments must be checked before trying to embed 
 *  structures with completely nonsensical constraints. Perhaps restrict A-A 
 *  ligands with bridge length 4 (chelating atoms included), maybe even up to 6
 *  to cis arrangements. Xantphos (with bridge length 7 is the smallest 
 *  trans-spanning ligand mentioned in Wikipedia).
 */
template<
  template<typename T = AssignmentColumn>
  class Symmetry
>
std::vector<
  Assignment<Symmetry>
> uniqueAssignments(
  const Assignment<Symmetry>& initial,
  const bool& removeTransSpanningGroups = true
) {
  std::vector<
    Assignment<Symmetry>
  > uniqueAssignments = {initial};

  auto rotationsSet = initial.generateAllRotations();

  Assignment<Symmetry> assignment = initial;

  do {
    bool currentAssignmentIsDistinct = rotationsSet.count(
      assignment
    ) == 0;
    if(currentAssignmentIsDistinct) {
      uniqueAssignments.push_back(assignment);
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
  } while(assignment.nextPermutation());

  if(removeTransSpanningGroups) {
    uniqueAssignments.erase(
      std::remove_if(
        uniqueAssignments.begin(),
        uniqueAssignments.end(),
        predicateHasTransArrangedPairs<Symmetry>
      ),
      uniqueAssignments.end()
    );
    return uniqueAssignments;
  } else return uniqueAssignments;
}

} // eo namespace

#endif
