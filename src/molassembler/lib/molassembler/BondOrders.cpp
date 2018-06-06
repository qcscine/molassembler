#include "BondOrders.h"
#include "BondDistance.h"

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
