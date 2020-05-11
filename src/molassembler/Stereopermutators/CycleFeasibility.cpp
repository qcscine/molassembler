/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Stereopermutators/CycleFeasibility.h"

#include "molassembler/Detail/CyclicPolygons.h"

namespace Scine {
namespace Molassembler {
namespace Stereopermutators {

namespace {

bool distanceToBaseContradictsGraph(
  const double cyclicPolygonDistance,
  const BaseAtom& base,
  const Utils::ElementType oppositeElement
) {
  // Add the out-of-plane distance to the cyclic polygon distance
  const double distance = std::sqrt(
    cyclicPolygonDistance * cyclicPolygonDistance
    + base.outOfPlaneDistance * base.outOfPlaneDistance
  );

  const double singleBondDistance = Bond::calculateBondDistance(
    base.elementType,
    oppositeElement,
    BondType::Single
  );

  // std::cout << distance << " <= " << singleBondDistance << " -> " << (distance <= singleBondDistance) << "\n";

  return distance <= singleBondDistance;
}

bool modelContradictsGraph(
  const std::vector<Utils::ElementType>& elementTypes,
  const std::vector<double>& cycleEdgeLengths,
  const std::vector<double>& phis,
  const BaseAtom& base
) {
  /* Layout of parameters:
   *
   * elementTypes: A, I, ..., X, B
   * cycleEdgeLengths: A-I, ..., X-B, B-A
   * phis: A-I-J, ..., X-B-A, B-A-I
   */
  const double alpha = CommonTrig::lawOfCosinesAngle(
    base.distanceToLeft,
    cycleEdgeLengths.back(), // B-A
    base.distanceToRight
  );

  const double d1 = CommonTrig::lawOfCosines(
    base.distanceToLeft, // Base-A
    cycleEdgeLengths.front(), // A-I
    alpha + phis.back() // Base-A-B + B-A-I
  );

  if(
    distanceToBaseContradictsGraph(
      d1,
      base,
      elementTypes.at(1) // I
    )
  ) {
    return true;
  }

  std::vector<double> distances {base.distanceToLeft, d1};
  std::vector<double> deltaAngles {};
  distances.reserve(phis.size());
  deltaAngles.reserve(phis.size());

  for(unsigned i = 0; i < phis.size() - 3; ++i) {
    /* Calculate a new delta angle so that we can construct the next triangle
     *
     * First edge: A-I, I-J, ...,
     * Second edge: d1, d2, ...,
     * Edge opposite to angle: Base-A, d1, ...,
     */
    deltaAngles.push_back(
      CommonTrig::lawOfCosinesAngle(
        cycleEdgeLengths.at(i),
        *(distances.end() - 1),
        *(distances.end() - 2)
      )
    );

    /* Calculate a new distance to the base
     *
     * First edge: d1, d2, ...,
     * Second edge: I-J, J-K, ...,
     * Angle: (A-I-J, I-J-K, ..., )
     *       - last delta angle
     */
    distances.push_back(
      CommonTrig::lawOfCosines(
        distances.back(),
        cycleEdgeLengths.at(i + 1),
        phis.at(i)
        - deltaAngles.back()
      )
    );

    // Check the distance against a hypothetical single bond distance
    if(
      distanceToBaseContradictsGraph(
        distances.back(),
        base,
        elementTypes.at(i + 2) // J, K, ...
      )
    ) {
      return true;
    }
  }

  return false;
}

} // namespace

bool cycleModelContradictsGraph(
  const std::vector<Utils::ElementType>& elementTypes,
  const std::vector<double>& cycleEdgeLengths,
  const std::vector<BaseAtom>& bases
) {
  if(!cyclic_polygons::exists(cycleEdgeLengths)) {
    return true;
  }

  assert(elementTypes.size() == cycleEdgeLengths.size());

  // std::cout << Temple::stringify(elementTypes) << " -> " << Temple::stringify(cycleEdgeLengths) << "\n";
  // const BaseAtom& frontBase = bases.front();
  // std::cout << "against distances " << frontBase.distanceToLeft << " and " << frontBase.distanceToRight << "\n";

  /* Calculate internal angles of the cyclic polygon. These are laid out:
   * elementTypes: A, I, ..., X, B
   * cycleEdgeLengths: A-I, ..., X-B, B-A
   * phis: A-I-J, ..., X-B-A, B-A-I
   */
  const auto phis = cyclic_polygons::internalAngles(cycleEdgeLengths);

  return Temple::any_of(
    bases,
    [&](const BaseAtom& base) -> bool {
      return modelContradictsGraph(elementTypes, cycleEdgeLengths, phis, base);
    }
  );
}

bool triangleBondTooClose(
  const double a,
  const double b,
  const double angle,
  const double bondRadius
) {
  // No need to calculate for angle ~= pi
  if(std::fabs(angle - M_PI) <= 1e-10) {
    return true;
  }

  /* We assume that the angles at the other points are acute, there would need
   * to be a pretty massive relative difference between a and b for one of them
   * to be obtuse and hence the altitude lying outside the triangle
   */
  const double c = CommonTrig::lawOfCosines(a, b, angle);
  const double gamma = CommonTrig::lawOfCosinesAngle(a, c, b);
  const double altitude = a * std::sin(gamma);

  return altitude <= bondRadius;
}

} // namespace Stereopermutators
} // namespace Molassembler
} // namespace Scine
