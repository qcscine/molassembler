#include "DelibHelpers.h"
#include "temple/Containers.h"

#include <Eigen/Geometry>

namespace molassembler {

namespace DelibHelpers {

bool validPositionIndices(
  const Delib::PositionCollection& positions,
  std::initializer_list<AtomIndexType> indices
) {
  /* scine::Delib casts its' underlying vector<Vector3>::size_type to int, so
   * we have to get at the unsigned std::size_t differently. Since vector
   * iterators are RandomAccessIterators, we can get a std::ptrdiff_t in O(1),
   * which is typically bigger than int in most data models.
   */
  auto signedPtrDiffDistance = std::distance(
    std::begin(positions),
    std::end(positions)
  );

  // Casting that to size_t gives us the maximum addressable space back
  std::size_t unsignedPositionCollectionSize = signedPtrDiffDistance;

  return temple::all_of(
    indices,
    [&unsignedPositionCollectionSize](const AtomIndexType& index) -> bool {
      return index < unsignedPositionCollectionSize;
    }
  );
}

double getDistance(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j
) {
  assert(
    validPositionIndices(positions, {i, j})
  );

  return (
    positions[i].asEigenVector()
    - positions[j].asEigenVector()
  ).norm();
}

double getAngle(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k
) {
  assert(
    validPositionIndices(positions, {i, j, k})
  );

  auto a = positions[i].asEigenVector() - positions[j].asEigenVector(),
       b = positions[k].asEigenVector() - positions[j].asEigenVector();

  return std::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

double getDihedral(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
) {
  assert(
    validPositionIndices(positions, {i, j, k, l})
  );

  Eigen::Vector3d a = positions[j].asEigenVector() - positions[i].asEigenVector();
  Eigen::Vector3d b = positions[k].asEigenVector() - positions[j].asEigenVector();
  Eigen::Vector3d c = positions[l].asEigenVector() - positions[k].asEigenVector();

  return std::atan2(
    (
      a.cross(b)
    ).cross(
      b.cross(c)
    ).dot(
      b.normalized()
    ),
    (
      a.cross(b)
    ).dot(
      b.cross(c)
    )
  );
}

double getDihedral(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndexType, 4>& indices
) {
  return getDihedral(
    positions,
    indices[0],
    indices[1],
    indices[2],
    indices[3]
  );
}


double getSignedVolume(
  const Delib::PositionCollection& positions,
  const AtomIndexType& i,
  const AtomIndexType& j,
  const AtomIndexType& k,
  const AtomIndexType& l
) {
  assert(
    validPositionIndices(positions, {i, j, k, l})
  );

  return (
    positions[i].asEigenVector()
    - positions[l].asEigenVector()
  ).dot(
    (
      positions[j].asEigenVector()
      - positions[l].asEigenVector()
    ).cross(
      positions[k].asEigenVector()
      - positions[l].asEigenVector()
    )
  );
}

double getSignedVolume(
  const Delib::PositionCollection& positions,
  const std::array<AtomIndexType, 4>& indices
) {
  return getSignedVolume(
    positions,
    indices[0],
    indices[1],
    indices[2],
    indices[3]
  );
}

} // namespace DelibHelpers

} // namespace molassembler
