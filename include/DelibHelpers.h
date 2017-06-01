#ifndef INCLUDE_DELIB_HELPERS_H
#define INCLUDE_DELIB_HELPERS_H

#include "Delib/PositionCollection.h"
#include "common_typedefs.h"

namespace MoleculeManip {

namespace DelibHelpers {

double getDistance(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j
);

double getAngle(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k
);

double getDihedral(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
);

double getDihedral(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndexType, 4>& indices
);


double getSignedVolume(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
);

double getSignedVolume(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndexType, 4>& indices
);

}

}

#endif
