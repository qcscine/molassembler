#include "BondDistance.h"

namespace molassembler {

namespace Bond {

double calculateBondDistance(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
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

} // namespace Bond

} // namespace molassembler
