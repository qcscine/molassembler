#include "chemical_symmetries/CoplanarityAnalysis.h"

#include <Eigen/SVD>

namespace Scine {
namespace Symmetry {

CoplanarityAnalysis analyseCoplanarGroups(Name symmetryName) {
  CoplanarityAnalysis analysis;

  analysis.coplanarGroups = detail::maximalCoplanarGroups(symmetryName);
  analysis.groupRelationships = detail::groupRelationships(
    symmetryName,
    analysis.coplanarGroups
  );

  return analysis;
}

namespace detail {

// TODO this does not include the midpoint at 0,0,0 yet at all, neither as sym pos nor as coordinate

Plane planeFit(
  const Name symmetryName,
  const std::vector<unsigned>& groupSymmetryPositions
) {
  const unsigned G = groupSymmetryPositions.size();
  assert(G > 3);

  const auto& symmetryCoordinates = symmetryData().at(symmetryName).coordinates;

  const Eigen::Vector3d centroid = temple::accumulate(
    groupSymmetryPositions,
    Eigen::Vector3d::Zero(),
    [&](const Eigen::Vector3d& carry, const unsigned symmetryPosition) -> Eigen::Vector3d {
      return carry + symmetryCoordinates.at(symmetryPosition);
    }
  ) / static_cast<double>(groupSymmetryPositions.size());

  const Eigen::Matrix<double, 3, Eigen::Dynamic> X (3, G);
  for(unsigned i = 0; i < G; ++i) {
    X.col(i) = symmetryCoordinates.at(
      groupSymmetryPositions.at(i)
    );
  }

  auto decomposition = X::jacobiSVD(Eigen::ComputeThinU | Eigen::ComputeThinV);

  // TODO continue
}

bool planeMembership(
  const Name symmetryName,
  const Plane& plane,
  const std::vector<unsigned>& symmetryPositions
) {
}

std::vector<Group> maximalCoplanarGroups(const Name symmetryName) {
  const unsigned S = size(symmetryName);

  for(unsigned groupSize = S; groupSize > 3; --groupSize) {
    std::vector<bool> mask (groupSize, true);
    mask.resize(S, false);

    do {
      std::vector<unsigned> pickedPositions;
      pickedPositions.reserve(groupSize);
      for(unsigned i = 0; i < S; ++i) {
        if(mask[i]) {
          pickedPositions.push_back(i);
        }
      }

      auto plane = planeFit(symmetryName, pickedPositions);

      // TODO All positions must fit closely to the plane
    } while(temple::inplace::next_permutation(mask));
  }
}

PlaneRelationship groupRelationship(
  Name symmetryName,
  const Group& a,
  const Group& b
) {
  const auto& coordinates = symmetryData().at(symmetryName).coordinates;

  /* Check the angle between normals to determine whether they are parallel,
   * perpendicular or unrelated
   */
  const auto& aNormal = a.plane.normal();
  const auto& bNormal = b.plane.normal();

  assert(aNormal.norm() == 1.0);
  assert(bNormal.norm() == 1.0);
  const double radianAngle = std::acos(
    aNormal.dot(bNormal)
  );

  constexpr double degreeTolerance = 2;

  if(
    radianAngle <= temple::Math::toRadians(degreeTolerance)
    || radianAngle >= temple::Math::toRadians(180 - degreeTolerance)
  ) {
    return PlaneRelationship::Parallel;
  }

  if(
    temple::Math::toRadians(90 - degreeTolerance) <= radianAngle
    && radianAngle <= temple::Math::toRadians(90 + degreeTolerance)
  ) {
    return PlaneRelationship::Perpendicular;
  }

  return PlaneRelationship::None;
}

std::vector<RelationshipTuple> groupRelationships(
  Name symmetryName,
  const std::vector<Group>& groups
) {
  std::vector<RelationshipTuple> relationships;

  for(unsigned i = 0; i < groups.size(); ++i) {
    for(unsigned j = i + 1; j < groups.size(); ++j) {
      PlaneRelationship relationship = groupRelationship(
        symmetryName,
        groups.at(i),
        groups.at(j)
      );

      if(relationship != PlaneRelationship::None) {
        relationships.push_back(i, j, relationship);
      }
    }
  }

  return relationships;
}

} // namespace detail

} // namespace Symmetry
} // namespace Scine
