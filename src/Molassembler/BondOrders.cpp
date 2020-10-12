/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/BondOrders.h"

#include "Utils/Typenames.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondDetector.h"

#include "Molassembler/AngstromPositions.h"
#include "Molassembler/Modeling/BondDistance.h"

namespace Scine {
namespace Molassembler {

Utils::BondOrderCollection uffBondOrders(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromPositions
) {
  const unsigned N = elements.size();

  Utils::BondOrderCollection bondOrders(N);

  for(unsigned i = 0; i < N; ++i) {
    for(unsigned j = i + 1; j < N; ++j) {
      double bondOrder = Bond::calculateBondOrder(
        elements.at(i),
        elements.at(j),
        (
          angstromPositions.positions.row(j)
          - angstromPositions.positions.row(i)
        ).norm()
      );

      if(bondOrder > 6.5) {
        throw std::logic_error(
          "Structure bond order interpretation yields bond orders greater than "
          "sextuple. The structure is most likely unreasonable."
        );
      }

      bondOrders.setOrder(i, j, bondOrder);
    }
  }

  return bondOrders;
}

Utils::BondOrderCollection covalentRadiiBondOrders(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromPositions
) {
  Utils::AtomCollection ac(elements, angstromPositions.getBohr());
  return Utils::BondDetector::detectBonds(ac);
}

} // namespace Molassembler
} // namespace Scine
