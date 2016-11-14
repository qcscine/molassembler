#ifndef INCLUDE_BOND_DISTANCE_H
#define INCLUDE_BOND_DISTANCE_H

#include <cmath>
#include <map>
#include "ElementTypes.h" // Delib

#include "AtomInfo.h"
#include "common_typedefs.h"
#include "Cache.h"

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

class DistanceCalculator {
private:
  using TupleType = std::tuple<
    Delib::ElementType, 
    Delib::ElementType, 
    BondType
  >;

  mutable std::map<
    TupleType,
    double
  > _calculated;

public:
  double get(
    const Delib::ElementType& a,
    const Delib::ElementType& b,
    const BondType& bondType
  ) const {
    auto tuple = TupleType(a, b, bondType);
    auto found = _calculated.find(tuple);

    if(found != _calculated.end()) {
      return found -> second;
    } else {
      double dist = calculateBondDistance(a, b, bondType);
      _calculated.insert(
        std::make_pair(
          tuple,
          dist
        )
      );
      return dist;
    }
  }
};

const DistanceCalculator distanceCalculator;

} // eo namespace Bond

} // eo namespace

#endif

