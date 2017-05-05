#include "BondDistance.h"

namespace MoleculeManip {

namespace Bond {

double calculateBondDistance(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) {
  return (
    AtomInfo::bondRadii.at(a)
    + AtomInfo::bondRadii.at(b)
    - ( // BO correction
      bondOrderCorrectionLambda * (
        AtomInfo::bondRadii.at(a)
        + AtomInfo::bondRadii.at(b)
      ) * log( bondOrderMap.at(static_cast<unsigned>(bondType)) )
    )
  );
}

} // namespace Bond

} // namespace MoleculeManip
