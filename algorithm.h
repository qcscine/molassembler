#ifndef LIB_UNIQUE_ASSIGNMENTS_ALGORITHM_HPP
#define LIB_UNIQUE_ASSIGNMENTS_ALGORITHM_HPP

#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
#include <iostream>
#include <sstream>

#include "Assignment.h"

/* TODO
 * - Bug in predicate* for monodentate ligands. For some reason, Assignments 
 *   are removed in post if the default on uniqueAssignments' removeTrans* is 
 *   true although no groups are present. Setting default to false leads to all
 *   monodentate tests passing.
 */

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
  const bool& removeTransSpanningGroups = false
) {
  std::vector<
    Assignment<Symmetry>
  > uniqueAssignments;

  Assignment<Symmetry> assignment = initial;

  do {
    bool currentAssignmentIsDistinct = std::accumulate(
      uniqueAssignments.begin(),
      uniqueAssignments.end(),
      true,
      [&assignment](
        const bool& carry,
        const Assignment<Symmetry>& uniqueAssignment
      ) {
        if(carry) { 
          // only bother to compute if the accumulation is still true
          return (
            carry
            && !uniqueAssignment.isRotationallySuperimposable(
              assignment
            )
          );
        } else {
          return carry;
        }
      }
    );
    if(currentAssignmentIsDistinct) {
      uniqueAssignments.push_back(assignment);
    }
  } while(assignment.nextPermutation());

  if(removeTransSpanningGroups) {
    uniqueAssignments.erase(
      std::remove_if(
        uniqueAssignments.begin(),
        uniqueAssignments.end(),
        predicateHasTransArrangedPairs<Symmetry>
      )
    );
    return uniqueAssignments;
  } else return uniqueAssignments;
}

#endif
