/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Modeling/BondDistance.h"

#include <cmath>

namespace Scine {

namespace molassembler {

namespace Bond {

double calculateBondDistance(
  const Utils::ElementType a,
  const Utils::ElementType b,
  const BondType bondType
) {
  return (
    atom_info::bondRadius(a)
    + atom_info::bondRadius(b)
    - ( // bond-order correction
      bondOrderCorrectionLambda * (
        atom_info::bondRadius(a)
        + atom_info::bondRadius(b)
      ) * log(
        bondOrderMap.at(
          static_cast<unsigned>(bondType)
        )
      )
    )
  );
}

double calculateBondOrder(
  const Utils::ElementType a,
  const Utils::ElementType b,
  const double distance
) {
  return std::exp(
    (
      atom_info::bondRadius(a) + atom_info::bondRadius(b) - distance
    ) / (
      bondOrderCorrectionLambda * (
        atom_info::bondRadius(a) + atom_info::bondRadius(b)
      )
    )
  );
}

} // namespace Bond

} // namespace molassembler

} // namespace Scine
