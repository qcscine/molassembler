#ifndef INCLUDE_DELIB_HELPERS_H
#define INCLUDE_DELIB_HELPERS_H

#include "Delib/PositionCollection.h"
#include "SharedTypes.h"

/*! @file
 *
 * A series of helpers to interface with the Delib library
 */

namespace molassembler {

namespace DelibHelpers {

//! Fetches the cartesian distance between two indices in a PositionCollection
double getDistance(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j
);

//! Fetches the angle in radians between three indices in a PostionCollection
double getAngle(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k
);

/*!
 * Fetches the dihedral angle in radians defined by four indices in a
 * PositionCollection.
 */
double getDihedral(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
);

/*!
 * Fetches the dihedral angle in radians defined by four indices in a
 * PositionCollection, but passed as an array
 */
double getDihedral(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndexType, 4>& indices
);

/*!
 * Returns the signed tetrahedron volume spanned by four indices in a
 * PositionCollection
 */
double getSignedVolume(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
);

/*!
 * Returns the signed tetrahedron volume spanned by four indices in a
 * PositionCollection, but passed as an array
 */
double getSignedVolume(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndexType, 4>& indices
);

/* Reimplementation on vector basis alone */

Eigen::Vector3d averagePosition(
  const Delib::PositionCollection& positions,
  const std::vector<AtomIndexType>& indices
);

double distance(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j
);

double angle(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k
);

double dihedral(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
);

/*!
 * Returns the signed tetrahedron volume spanned by four spatial positions
 * adjusted by V' = 6 * V
 */
double adjustedSignedVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
);

}

}

#endif
