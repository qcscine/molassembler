/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_CYCLE_FEASIBILITY_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATORS_CYCLE_FEASIBILITY_H

#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/Types.h"

#include "CyclicPolygons.h"

namespace Scine {
namespace molassembler {
namespace Stereopermutators {

struct BaseAtom {
  Utils::ElementType elementType;
  double distanceToLeft = 0;
  double distanceToRight = 0;
  double outOfPlaneDistance = 0;
};

/**
 * @brief Checks whether spatially modeling a cycle contradicts the graph
 *
 * @param elementTypes Element types of the cycle, laid out as A, I, J, ..., X, B
 * @param cycleEdgeLengths Edge lengths forming a cyclic polygon, laid out as
 *   A-I, I-J, ..., X-B, B-A
 * @param bases Base atoms to check distances to atoms of the cyclic polygon
 *   against. Each base is distanceToLeft away from A and distanceToRight away
 *   from B
 *
 * @return If any distances to atoms of the cyclic polygon are shorter
 *   than a single bond
 */
bool cycleModelContradictsGraph(
  const std::vector<Utils::ElementType>& elementTypes,
  const std::vector<double>& cycleEdgeLengths,
  const std::vector<BaseAtom>& bases
);

/**
 * @brief Checks whether the far bond in a triangle is too close to the
 *   originating atom
 *
 * @param a Bond length to first atom of triangle
 * @param b Bond length to second atom of triangle
 * @param angle Angle between bonds @p a and @p b
 * @param bondRadius Bond radius of the originating atom
 *
 * @returns Whether the altitude of the originating atom in the resulting
 *   triangle is shorter or equal to the bond radius
 */
bool triangleBondTooClose(
  const double a,
  const double b,
  const double angle,
  const double bondRadius
);

} // namespace Stereopermutators
} // namespace molassembler
} // namespace Scine

#endif
