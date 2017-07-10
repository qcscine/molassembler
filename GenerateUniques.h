#ifndef LIB_UNIQUE_ASSIGNMENTS_GENERATE_UNIQUES_H
#define LIB_UNIQUE_ASSIGNMENTS_GENERATE_UNIQUES_H

#include <vector>

#include "Assignment.h"

/*! @file
 *
 * Main entry point into the library. From here, a set of rotationally unique
 * assignments for a specific assignment can be generated, with or without
 * counting the relative weights of unique assignments.
 */

namespace UniqueAssignments {

bool predicateHasTransArrangedPairs(
  const Assignment& assignment
);

/*! 
 * Generate the set of rotationally unique assignments for a given assignment.
 * By default removes trans-spanning groups (where a linked group's
 * directly bonded atoms span an angle of 180°).
 *
 * NOTE: Gives NO guarantees as to satisfiability (if assignments can be
 *  fulfilled with real ligands) 
 * E.g. M (A-A)_3 generates a trans-trans-trans assignment, which is extremely 
 *  hard to find actual ligands for that work.
 * The satisfiability of assignments must be checked before trying to embed 
 *  structures with completely nonsensical constraints. Perhaps restrict A-A 
 *  ligands with bridge length 4 (chelating atoms included), maybe even up to 6
 *  to cis arrangements. Xantphos (with bridge length 7) is the smallest 
 *  trans-spanning ligand mentioned in Wikipedia.
 */
std::vector<Assignment> uniqueAssignments(
  const Assignment& initial,
  const bool& removeTransSpanningGroups = true
);

//! Data class for uniqueAssignments including weights
struct UniqueAssignmentsReturnType {
  std::vector<Assignment> assignments;
  std::vector<unsigned> weights;
};

/*! 
 * Returns the set of rotationally unique assignments including absolute
 * occurrence counts.
 */
UniqueAssignmentsReturnType uniqueAssignmentsWithCounts(
  const Assignment& initial,
  const bool& removeTransSpanningGroups = true
);

} // eo namespace

#endif
