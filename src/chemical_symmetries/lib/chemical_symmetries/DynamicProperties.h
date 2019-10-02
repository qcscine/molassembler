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

namespace Scine {

namespace Symmetry {

namespace properties {

constexpr double floatingPointEqualityThreshold [[gnu::unused]] = 1e-4;

/*! @brief Rotates a passed list of indices with a specified rotation vector
 *
 * @complexity{@math{\Theta(S)}}
 */
std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const std::vector<unsigned>& rotation
);

/*! @brief Rotates a passed list of indices of a specific symmetry
 *
 * @complexity{@math{\Theta(S)}}
 */
std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  Name symmetryName,
  unsigned rotationFunctionIndex
);

/*! @brief Calculate the periodicty of a symmetry's rotation
 *
 * @complexity{@math{\Theta(M S)} where @math{M} is the multiplicity of the
 * rotation and @math{S} is the symmetry size}
 */
unsigned rotationPeriodicity(
  Name symmetryName,
  const std::vector<unsigned>& rotation
);

/*! @brief Generate a character representation of a symmetry's position groups
 *
 * Groups symmetry positions according to matching sets of angles towards all
 * other symmetry positions. Then transforms the symmetry positions themselves
 * to character representations of the groups they belong to.
 *
 * @complexity{@math{\Theta(S^2)}}
 *
 * @param symmetryName name of the symmetry for which to generate the character
 *   representation
 */
std::vector<char> positionGroups(Name symmetryName);

/*! @brief Generate the inverse rotation to a symmetry's rotation
 *
 * Inverts the permutation provided by @p rotation. UB if rotation is not an
 * index permutation.
 *
 * @complexity{@math{\Theta(N)}}
 */
std::vector<unsigned> inverseRotation(const std::vector<unsigned>& rotation);

/*! @brief Gets the coordinates of an indexOptional for a specific symmetry.
 *
 * As was defined previously, boost::none is a placeholder for the central atom,
 * which is not explicitly held in memory as it is always placed at {0, 0, 0}.
 *
 * @complexity{@math{\Theta(1)}}
 */
Eigen::Vector3d getCoordinates(
  Name symmetryName,
  const boost::optional<unsigned>& indexInSymmetryOption
);

/*! @brief Tetrahedron volume spanned by four positions
 *
 * Returns the Eigen calculation of a signed tetrahedron volume from four
 * tetrahedron edge point vectos
 *
 * @complexity{@math{\Theta(1)}}
 */
double getTetrahedronVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
);

/*! @brief Calculates angular distorition for an index mapping between symmetries
 *
 * Calculates the angular distortion for a transition between two symmetries and
 * a specific index mapping that connects indices from the source symmetry to
 * the target symmetry
 *
 * @complexity{@math{\Theta(S^2)}}
 */
double calculateAngleDistortion(
  Name from,
  Name to,
  const std::vector<unsigned>& indexMapping
);

/*! @brief Propagates an index optional through an index mapping
 *
 * Index optionals as used in tetrahedron definitions in symmetries need to be
 * propagated through index mappings generated in the process of finding the
 * optimal symmetry transitions. Shorthand function that performs a propagation
 * if a the passed indexOptional is not boost::none, returns none otherwise.
 *
 * @complexity{@math{\Theta(1)}}
 */
boost::optional<unsigned> propagateIndexOptionalThroughMapping(
  const boost::optional<unsigned>& indexOptional,
  const std::vector<unsigned>& indexMapping
);

/*! @brief Calculates chiral distortion for a transition between symmetries
 *
 * Calculates the chiral distortion for a transition between two symmetries and
 * a specific index mapping that connects indices from the source symmetry to
 * the target symmetry
 *
 * @complexity{@math{\Theta(T)} where @math{T} is the number of tetrahedra for
 * the symmetry, typically small}
 */
double calculateChiralDistortion(
  Name from,
  Name to,
  const std::vector<unsigned>& indexMapping
);

/*! @brief Generates all rotations of a sequence of indices within a symmetry group
 *
 * @complexity{At most maxRotation iterations}
 */
std::set<
  std::vector<unsigned>
> generateAllRotations(
  Name symmetryName,
  const std::vector<unsigned>& indices
);

/*! @brief Transform symmetry positions through a mapping
 *
 * Writes the indices of the original symmetry in the mapping into the target
 * symmetry's indexing scheme.
 *
 * @complexity{@math{\Theta(S)}}
 *
 * @param mapping An index mapping that specifies how indices are mapped
 *   from a source symmetry to a target symmetry
 */
std::vector<unsigned> applyIndexMapping(
  Name to,
  const std::vector<unsigned>& mapping
);

/**
 * @brief Data type grouping distortions between symmetries
 */
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

/*! @brief Calculates ideal index mappings for +1, 0 size transitions
 *
 * Generates symmetry transition index mappings with the lowest angular
 * distortion and then subsets that group to those with the lowest chiral
 * distortion. Transitions are limited to symmetries with size differences of 0
 * and ±1.
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @param symmetryFrom Transition source symmetry
 * @param symmetryTo Transition target symmetry
 *
 * @pre Symmetry::size(symmetryFrom) + {0, 1} == Symmetry::size(symmetryTo)
 */
std::vector<DistortionInfo> symmetryTransitionMappings(
  Name symmetryFrom,
  Name symmetryTo
);

/*! @brief Calculates ideal index mappings for ligand loss transitions
 *
 * Generates symmetry transition index mappings for the special case of ligand
 * loss, in which a ligand is removed from a particular position in the symmetry
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @param symmetryFrom Transition source symmetry
 * @param symmetryTo Transition target symmetry
 * @param positionInSourceSymmetry position lost in source symmetry
 *
 * @pre Symmetry::size(symmetryFrom) == Symmetry::size(symmetryTo) + 1
 */
std::vector<DistortionInfo> ligandLossTransitionMappings(
  Name symmetryFrom,
  Name symmetryTo,
  unsigned positionInSourceSymmetry
);

//! A grouping of index mappings of equal angular and chiral distortion
struct SymmetryTransitionGroup {
  /*!
   * @brief A list of index mappings that share the same @p angularDistortion
   *   and @p chiralDistortion
   */
  std::vector<
    std::vector<unsigned>
  > indexMappings;
  double angularDistortion;
  double chiralDistortion;

  SymmetryTransitionGroup(
    std::vector<
      std::vector<unsigned>
    > passIndexMappings,
    const double& passAngleDistortion,
    const double& passChiralDistortion
  );

  SymmetryTransitionGroup() = default;
};

/*! @brief Selects the best transition mapping from many DistortionInfos
 *
 * Selects index mappings from a DistortionInfo list, choosing those with lowest
 * angular distortion first, and lowest chiral distortion afterwards. Groups
 * them into a SymmetryTransitionGroup.
 *
 * @complexity{@math{\Theta(N)}}
 */
SymmetryTransitionGroup selectBestTransitionMappings(
  const std::vector<DistortionInfo>& distortions
);

/*!
 * @brief Calculates the number of stereopermutations in a specific symmetry
 *   and a number of identical ligands.
 *
 * @complexity{@math{\Theta(S!)}}
 */
unsigned numUnlinkedStereopermutations(
  Name symmetry,
  unsigned nIdenticalLigands
);

/*!
 * @brief Calculates if there are multiple unlinked stereopermutations in a
 *   specific symmetry for a number of identical ligands.
 *
 * @complexity{@math{\Theta(S!)}}
 */
bool hasMultipleUnlinkedStereopermutations(
  Name symmetry,
  unsigned nIdenticalLigands
);

/*! @brief Yields the symmetry with the most rotations from a selection
 *
 * @complexity{@math{\Theta(1)}}
 */
Name mostSymmetric(std::vector<Name> selection);

/*! @brief Yields the symmetry with the most rotations of a particular size
 *
 * @complexity{@math{\Theta(1)}}
 */
Name mostSymmetric(unsigned symmetrySize);

} // namespace properties

} // namespace Symmetry

} // namespace Scine

#endif
