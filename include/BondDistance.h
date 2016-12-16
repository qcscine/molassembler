#ifndef INCLUDE_BOND_DISTANCE_H
#define INCLUDE_BOND_DISTANCE_H

#include <cmath>
#include <map>
#include "ElementTypes.h" // Delib

#include "AtomInfo.h"
#include "common_typedefs.h"

namespace MoleculeManip {

namespace Bond {

extern const std::map<BondType, double> bondOrderMap;

const double bondOrderCorrectionLambda = 0.1332;

double calculateBondDistance(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
);

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

extern const DistanceCalculator distanceCalculator;

} // eo namespace Bond

} // eo namespace

#endif

