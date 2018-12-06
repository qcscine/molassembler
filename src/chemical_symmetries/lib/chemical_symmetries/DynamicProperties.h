/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Run-time symmetry property calculations
 *
 * Contains a suite of property calculations on the dynamic symmetry data.
 */

#ifndef INCLUDE_SYMMETRY_DYNAMIC_PROPERTIES_CALCULATION_H
#define INCLUDE_SYMMETRY_DYNAMIC_PROPERTIES_CALCULATION_H

#include "Eigen/Core"
#include "boost/optional/optional_fwd.hpp"

#include "chemical_symmetries/Names.h"

#include <set>
#include <vector>


namespace Symmetry {

namespace properties {

constexpr double floatingPointEqualityThreshold [[gnu::unused]] = 1e-4;

//! Rotates a passed list of indices with a specified rotation vector
std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const std::vector<unsigned>& rotation
);

//! Rotates a passed list of indices of a specific symmetry
std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const Symmetry::Name symmetryName,
  unsigned rotationFunctionIndex
);

//! Calculate the periodicty of a symmetry's rotation
unsigned rotationPeriodicity(
  const Symmetry::Name symmetryName,
  const std::vector<unsigned>& rotation
);

//! Generate a character representation of a symmetry's position groups
std::vector<char> positionGroups(const Symmetry::Name symmetryName);

//! Generate the inverse rotation to a symmetry's rotation
std::vector<unsigned> inverseRotation(const std::vector<unsigned>& rotation);

//! Generate a symmetry's occupation sequence that causes all tetrahedra to invert
std::vector<unsigned> invertedSequence(const Symmetry::Name symmetryName);

/*!
 * Gets the coordinates of an indexOptional for a specific symmetry.
 * As was defined previously, boost::none is a placeholder for the central atom,
 * which is not explicitly held in memory as it is always placed at {0, 0, 0}.
 */
Eigen::Vector3d getCoordinates(
  const Symmetry::Name symmetryName,
  const boost::optional<unsigned>& indexInSymmetryOption
);

/*!
 * Returns the Eigen calculation of a signed tetrahedron volume from four
 * tetrahedron edge point vectos
 */
double getTetrahedronVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
);

/*!
 * Calculates the angular distortion for a transition between two symmetries and
 * a specific index mapping that connects indices from the source symmetry to
 * the target symmetry
 */
double calculateAngleDistortion(
  const Symmetry::Name from,
  const Symmetry::Name to,
  const std::vector<unsigned>& indexMapping
);

/*!
 * Index optionals as used in tetrahedron definitions in symmetries need to be
 * propagated through index mappings generated in the process of finding the
 * optimal symmetry transitions. Shorthand function that performs a propagation
 * if a the passed indexOptional is not boost::none, returns that otherwise.
 */
boost::optional<unsigned> propagateIndexOptionalThroughMapping(
  const boost::optional<unsigned>& indexOptional,
  const std::vector<unsigned>& indexMapping
);

/*!
 * Calculates the chiral distortion for a transition between two symmetries and
 * a specific index mapping that connects indices from the source symmetry to
 * the target symmetry
 */
double calculateChiralDistortion(
  const Symmetry::Name from,
  const Symmetry::Name to,
  const std::vector<unsigned>& indexMapping
);

/*!
 * Generates a set of all possible rotations of an index sequence within a
 * desired symmetry
 */
std::set<
  std::vector<unsigned>
> generateAllRotations(
  const Symmetry::Name symmetryName,
  const std::vector<unsigned>& indices
);

/*!
 * Writes the indices of the original symmetry in the mapping into the target
 * symmetry's indexing scheme.
 */
std::vector<unsigned> applyIndexMapping(
  const Symmetry::Name to,
  const std::vector<unsigned>& mapping
);

struct DistortionInfo {
  std::vector<unsigned> indexMapping;
  double angularDistortion;
  double chiralDistortion;

  DistortionInfo(
    std::vector<unsigned> passIndexMapping,
    const double& passAngularDistortion,
    const double& passChiralDistortion
  );
};

/*!
 * Generates symmetry transition index mappings with the lowest angular
 * distortion and then subsets that group to those with the lowest chiral
 * distortion. Transitions are limited to symmetries with size differences of 0
 * and ±1.
 */
std::vector<DistortionInfo> symmetryTransitionMappings(
  const Symmetry::Name symmetryFrom,
  const Symmetry::Name symmetryTo
);

/*!
 * Generates symmetry transition index mappings for the special case of ligand
 * loss, in which a ligand is removed from a particular position in the symmetry
 */
std::vector<DistortionInfo> ligandLossTransitionMappings(
  const Symmetry::Name symmetryFrom,
  const Symmetry::Name symmetryTo,
  const unsigned positionInSourceSymmetry
);

//! A grouping of index mappings of equal angular and chiral distortion
struct SymmetryTransitionGroup {
  std::vector<
    std::vector<unsigned>
  > indexMappings;
  double angularDistortion, chiralDistortion;

  SymmetryTransitionGroup(
    std::vector<
      std::vector<unsigned>
    > passIndexMappings,
    const double& passAngleDistortion,
    const double& passChiralDistortion
  );

  SymmetryTransitionGroup() {}
  SymmetryTransitionGroup(SymmetryTransitionGroup&& other)
    : indexMappings(std::move(other.indexMappings)),
      angularDistortion(other.angularDistortion),
      chiralDistortion(other.chiralDistortion)
  {}

  SymmetryTransitionGroup(const SymmetryTransitionGroup& other)
    : indexMappings(other.indexMappings),
      angularDistortion(other.angularDistortion),
      chiralDistortion(other.chiralDistortion)
  {}

  SymmetryTransitionGroup& operator = (const SymmetryTransitionGroup& other) {
    indexMappings = other.indexMappings;
    angularDistortion = other.angularDistortion;
    chiralDistortion = other.chiralDistortion;

    return *this;
  }
};

/*!
 * Selects index mappings from a DistortionInfo list, choosing those with lowest
 * angular distortion first, and lowest chiral distortion afterwards. Groups
 * them into a SymmetryTransitionGroup.
 */
SymmetryTransitionGroup selectBestTransitionMappings(
  const std::vector<DistortionInfo>& distortions
);

/*!
 * Calculates the number of assignments in a specific symmetry and a number of
 * identical ligands.
 */
unsigned numUnlinkedAssignments(
  const Symmetry::Name symmetry,
  const unsigned nIdenticalLigands
);

/*!
 * Calculates if there are multiple unlinked assignments in a specific symmetry
 * for a number of identical ligands.
 */
bool hasMultipleUnlinkedAssignments(
  const Symmetry::Name symmetry,
  const unsigned nIdenticalLigands
);

} // namespace properties

} // namespace Symmetry

#endif
