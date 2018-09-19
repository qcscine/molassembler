// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/BondOrders.h"

#include "Delib/ElementTypeCollection.h"
#include "Delib/BondOrderCollection.h"

#include "molassembler/AngstromWrapper.h"
#include "molassembler/Modeling/BondDistance.h"

namespace molassembler {

Delib::BondOrderCollection uffBondOrders(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper
) {
  const int N = elements.size();

  auto bondOrders = Delib::BondOrderCollection::createEmpty(N);

  // Integers to comply with Delib
  for(int i = 0; i < N; ++i) {
    for(int j = i + 1; j < N; ++j) {
      bondOrders.setOrder(
        i,
        j,
        Bond::calculateBondOrder(
          elements.at(i),
          elements.at(j),
          (
            angstromWrapper.positions.at(j)
            - angstromWrapper.positions.at(i)
          ).norm()
        )
      );
    }
  }

  return bondOrders;
}

} // namespace molassembler
