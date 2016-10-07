#ifndef INCLUDE_COMMON_TRIG_H
#define INCLUDE_COMMON_TRIG_H

#include <cmath>
#include "ElementTypes.h" // Delib

#include "BondDistance.h"

namespace MoleculeManip {

namespace CommonTrig {

template<typename T>
T lawOfCosines(
  const T& a,
  const T& b,
  const T& phi_degrees
) {
  return sqrt(
    a * a
    + b * b
    - 2 * a * b * cos( // convert to radians
      M_PI * phi_degrees / ( 
        180.0
      )
    )
  );
}

std::pair<double, double> oneTwoDistanceBounds(
  const std::pair<Delib::ElementType, Delib::ElementType>& elementTypes,
  const double& distanceVariation,
  const BondType& bondType
) {
  double bondDistance = Bond::calculateBondDistance(
    elementTypes.first,
    elementTypes.second,
    bondType
  );
  return std::make_pair(
    bondDistance - distanceVariation,
    bondDistance + distanceVariation
  );
}



} // eo CommonTrig

}

#endif
