/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Analyze symmetries for coplanarity of symmetry positions and their
 *   mutual relationships
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_COPLANARITY_ANALYSIS_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_COPLANARITY_ANALYSIS_H

#include "chemical_symmetries/Names.h"
#include <Eigen/Geometry>
#include <vector>
#include <tuple>

namespace Scine {
namespace Symmetry {

struct CoplanarityAnalysis {
  //! Plane representation
  using Plane = Eigen::Hyperplane<double, 3>;

  //! A group is an unordered list of symmetry positions
  struct Group {
    std::vector<unsigned> symmetryPositions;
    Plane plane;
  };

  /**
   * @brief Types of geometric relationships between coplanar groups
   */
  enum class PlaneRelationship {
    Perpendicular,
    Parallel,
    None
  };

  //! Relationships between coplanar groups via two group indices and a relationship enum
  using RelationshipTuple = std::tuple<unsigned, unsigned, PlaneRelationship>;

  //! A list of coplanar symmetry position groups
  std::vector<Group> coplanarGroups;

  //! A list of relationships between the coplanar groups.
  std::vector<RelationshipTuple> groupRelationships;
};

CoplanarityAnalysis analyseCoplanarGroups(Name symmetryName);

namespace detail {

/**
 * @brief Determine the groups of larges size that are approximately coplanar
 */
std::vector<Group> maximalCoplanarGroups(Name symmetryName);

PlaneRelationship groupRelationship(
  Name symmetryName,
  const Group& a,
  const Group& b
);

std::vector<RelationshipTuple> groupRelationships(
  Name symmetryName,
  const std::vector<Group>& groups
);

} // namespace detail

} // namespace Symmetry
} // namespace Scine

#endif
