#ifndef LIB_UNIQUE_ASSIGNMENTS_GENERATE_UNIQUES_H
#define LIB_UNIQUE_ASSIGNMENTS_GENERATE_UNIQUES_H

#include <vector>

#include "Assignment.h"

/* TODO
 */

namespace UniqueAssignments {

bool predicateHasTransArrangedPairs(
  const Assignment& assignment
);

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
std::vector<
  Assignment
> uniqueAssignments(
  const Assignment& initial,
  const bool& removeTransSpanningGroups = true
);

} // eo namespace

#endif
