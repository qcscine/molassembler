#include "BondDistance.h"

namespace MoleculeManip {

namespace Bond {

const DistanceCalculator distanceCalculator;

const std::map<BondType, double> bondOrderMap({
  {BondType::Single, 1},
  {BondType::Double, 2},
  {BondType::Triple, 3},
  {BondType::Quadruple, 4},
  {BondType::Quintuple, 5},
  {BondType::Sextuple, 6},
  {BondType::Aromatic,  1.5}
});

double calculateBondDistance(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) {
  return (
    AtomInfo::bondRadii.at(a)
    + AtomInfo::bondRadii.at(b)
    - ( // BO correction
      bondOrderCorrectionLambda
      * (
        AtomInfo::bondRadii.at(a)
        + AtomInfo::bondRadii.at(b)
      )
      * log( bondOrderMap.at(bondType) )
    )
  );
}

} // eo namespace Bond

} // eo namespace
