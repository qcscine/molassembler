#ifndef INCLUDE_BOND_DISTANCE_H
#define INCLUDE_BOND_DISTANCE_H

#include <cmath>
#include "ElementTypes.h" // Delib

#include "AtomInfo.h"
#include "common_typedefs.h"

namespace MoleculeManip {

namespace Bond {

const std::map<BondType, double> bondOrderMap({
  {BondType::Single, 1},
  {BondType::Double, 2},
  {BondType::Triple, 3},
  {BondType::Quadruple, 4},
  {BondType::Quintuple, 5},
  {BondType::Sextuple, 6},
  {BondType::Aromatic,  1.5}
});

const double bondOrderCorrectionLambda = 0.1332;

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

#endif

