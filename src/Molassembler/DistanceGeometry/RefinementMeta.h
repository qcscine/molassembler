/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Contains functions overarching different refinement implementations
 */

#ifndef INCLUDE_MOLASSEMBLER_REFINEMENT_META_H
#define INCLUDE_MOLASSEMBLER_REFINEMENT_META_H

#include "Molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "Molassembler/Log.h"

#include "Molassembler/Temple/Stringify.h"

namespace Scine {
namespace Molassembler {
namespace DistanceGeometry {

/**
 * @brief Decides whether the final structure from a refinement is acceptable
 *
 * A final structure is acceptable if
 * - All distance bounds are within 0.5 of either the lower or upper boundary
 * - All chiral constraints are within 0.5 of either the lower or upper boundary
 *
 * @complexity{@math{\Theta(N^2)} due to atom-pairwise distance bounds}
 *
 * @param refinement The refinement functor
 * @param bounds The distance bounds
 * @param positions The final positions from a refinement
 *
 * @return Whether the final structure is acceptable
 */
template<class RefinementType, typename PositionType>
bool finalStructureAcceptable(
  const RefinementType& refinement,
  const DistanceBoundsMatrix& bounds,
  const PositionType& positions
) {
  struct FinalStructureAcceptableVisitor {
    const double deviationThreshold = 0.5;
    bool earlyExit = false;
    bool value = true;

    void distanceOverThreshold(AtomIndex /* i */, AtomIndex /* j */, double /* distance */) {
      earlyExit = true;
      value = false;
    }

    void chiralOverThreshold(const ChiralConstraint& chiral, double /* volume */) {
      if(!chiral.targetVolumeIsZero()) {
        earlyExit = true;
        value = false;
      }
    }

    void dihedralOverThreshold(const DihedralConstraint& /* dihedral */, double /* angle */, double /* term */) {
      earlyExit = true;
      value = false;
    }
  };

  return refinement.visitUnfulfilledConstraints(
    bounds,
    positions,
    FinalStructureAcceptableVisitor {}
  );
}

} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine

#endif
