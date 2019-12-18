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
    AtomInfo::bondRadius(a)
    + AtomInfo::bondRadius(b)
    - ( // bond-order correction
      bondOrderCorrectionLambda * (
        AtomInfo::bondRadius(a)
        + AtomInfo::bondRadius(b)
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
      AtomInfo::bondRadius(a) + AtomInfo::bondRadius(b) - distance
    ) / (
      bondOrderCorrectionLambda * (
        AtomInfo::bondRadius(a) + AtomInfo::bondRadius(b)
      )
    )
  );
}

} // namespace Bond

} // namespace molassembler

} // namespace Scine
