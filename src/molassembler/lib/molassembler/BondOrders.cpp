// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/BondOrders.h"

#include "Utils/Typenames.h"
#include "Utils/Bonds/BondOrderCollection.h"

#include "molassembler/AngstromWrapper.h"
#include "molassembler/Modeling/BondDistance.h"

namespace molassembler {

Scine::Utils::BondOrderCollection uffBondOrders(
  const Scine::Utils::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper
) {
  const unsigned N = elements.size();

  Scine::Utils::BondOrderCollection bondOrders(N);

  for(unsigned i = 0; i < N; ++i) {
    for(unsigned j = i + 1; j < N; ++j) {
      double bondOrder = Bond::calculateBondOrder(
        elements.at(i),
        elements.at(j),
        (
          angstromWrapper.positions.row(j)
          - angstromWrapper.positions.row(i)
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

} // namespace molassembler
