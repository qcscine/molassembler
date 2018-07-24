#include "Composites.h"

#include "Eigen/Geometry"

#include "chemical_symmetries/DynamicProperties.h"
#include "temple/Containers.h"
#include "temple/Stringify.h"

#include <iostream>

namespace stereopermutation {

namespace detail {

void rotateCoordinates(
  std::vector<Eigen::Vector3d>& positions,
  const Eigen::Vector3d& unitSource,
  const Eigen::Vector3d& unitTarget
) {
  // Special case: unitTarget = unitSource. Done.
  if(unitSource == unitTarget) {
    return;
  }

  /* Special case: unitTarget = -unitSource. Here, the rotation is not uniquely
   * defined, so we just invert all positions instead.
   */
  if(unitSource == -unitTarget) {
    for(auto& position : positions) {
      position *= -1;
    }

    return;
  }

  // Adapted from https://math.stackexchange.com/q/476311
  // Cross product of the unit vectors
  Eigen::Vector3d v = unitSource.cross(unitTarget);

  // Skew-symmetric cross product matrix
  Eigen::Matrix3d v_x;
  v_x <<     0, -v.z(),  v.y(),
         v.z(),      0, -v.x(),
        -v.y(),  v.x(),      0;

  // Dot product (essentially the cosine of the angle for these unit vectors)
  double c = unitSource.dot(unitTarget);

  // Calculate the rotation matrix
  Eigen::Matrix3d rotation;
  rotation = Eigen::Matrix3d::Identity() + v_x + v_x * v_x * (1.0 / (1 + c));

  for(auto& position : positions) {
    position = rotation * position;
  }
}

void translateCoordinates(
  std::vector<Eigen::Vector3d>& positions,
  const Eigen::Vector3d& translation
) {
  for(auto& position : positions) {
    position += translation;
  }
}

/*! Calculates the dihedral between four positions
 *
 * \note Resulting dihedrals are distributed on (-M_PI, M_PI].
 */
double dihedral(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
) {
  Eigen::Vector3d a = j - i;
  Eigen::Vector3d b = k - j;
  Eigen::Vector3d c = l - k;

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

} // namespace detail

constexpr temple::floating::ExpandedAbsoluteEqualityComparator<double> Composite::fpComparator;

double Composite::perpendicularSubstituentAngle(
  const double angleFromBoundSymmetryPosition,
  const double angleBetweenSubstituents
) {
  assert(angleFromBoundSymmetryPosition != M_PI);

  return std::acos(
    1.0 - (
      1.0 - std::cos(angleBetweenSubstituents)
    ) / (
      std::pow(
        std::sin(angleFromBoundSymmetryPosition),
        2
      )
    )
  );
}

std::vector<unsigned> Composite::identityRotation(const Symmetry::Name symmetryName) {
  return temple::iota<unsigned>(Symmetry::size(symmetryName));
}

std::vector<unsigned> Composite::generateRotation(
  const Symmetry::Name symmetryName,
  const unsigned fixedSymmetryPosition,
  const std::vector<unsigned>& changedPositions
) {
  auto periodicities = temple::map(
    temple::iota<unsigned>(Symmetry::rotations(symmetryName).size()),
    [&symmetryName](const unsigned rotationFunctionIndex) -> unsigned {
      return Symmetry::properties::rotationPeriodicity(
        symmetryName,
        Symmetry::rotations(symmetryName).at(rotationFunctionIndex)
      );
    }
  );

  auto rotationAltersPositions = [&](const std::vector<unsigned>& rotation) -> bool {
    return temple::all_of(
      changedPositions,
      [&rotation](const unsigned symmetryPosition) -> bool {
        return rotation.at(symmetryPosition) != symmetryPosition;
      }
    );
  };

  // Which rotation indicates a rotation around the bound symmetry position?
  std::vector<unsigned> rotationUses (periodicities.size(), 0);
  ++rotationUses.back();

  std::vector<unsigned> rotation;
  bool rotationFound = false;

  do {
    std::vector<unsigned> rotationIndexApplicationSequence;
    rotationIndexApplicationSequence.reserve(30);
    for(unsigned r = 0; r < rotationUses.size(); ++r) {
      // Add r as many times as rotationUses says to the application sequence
      for(unsigned repeat = 0; repeat < rotationUses.at(r); ++repeat) {
        rotationIndexApplicationSequence.push_back(r);
      }
    }

    temple::sort(rotationIndexApplicationSequence);

    do {
      // Create the rotation using the index application sequence front-to-back
      rotation = identityRotation(symmetryName);

      for(const auto r : rotationIndexApplicationSequence) {
        rotation = Symmetry::properties::applyRotation(
          rotation,
          symmetryName,
          r
        );
      }

      if(
        rotation.at(fixedSymmetryPosition) == fixedSymmetryPosition
        && rotationAltersPositions(rotation)
      ) {
        rotationFound = true;
      }
    } while(
      !rotationFound
      && std::next_permutation(
        std::begin(rotationIndexApplicationSequence),
        std::end(rotationIndexApplicationSequence)
      )
    );

  } while(
    !rotationFound
    && detail::nextCombinationPermutation(rotationUses, periodicities)
  );

  if(rotationFound) {
    return rotation;
  }

  return {};
}

std::vector<unsigned> Composite::rotation(
  const Symmetry::Name symmetryName,
  const unsigned fixedSymmetryPosition,
  const std::vector<unsigned>& perpendicularPlanePositions
) {
  /* Three possibilities:
   */

  if(perpendicularPlanePositions.size() > 1) {
    /* There are multiple elements in perpendicularPlanePositions. We have to
     * generate a rotation that keeps fixedSymmetryPosition fixed but rotates the
     * perpendicularPlanePositions, ideally with a periodicity equivalent to the
     * amount of symmetry positions involved.
     */
    auto candidateRotation = generateRotation(
      symmetryName,
      fixedSymmetryPosition,
      perpendicularPlanePositions
    );

    // There may be multiple elements, but no rotation. Return identity
    if(candidateRotation.empty()) {
      return {1};
    }

    /* Require that the periodicity of the discovered rotation is equal to the
     * number of elements being rotated. This should be a natural property of
     * generateRotation, but it's best to be sure.
     */
    assert(
      Symmetry::properties::rotationPeriodicity(symmetryName, candidateRotation)
      == perpendicularPlanePositions.size()
    );

    return candidateRotation;
  }

  if(perpendicularPlanePositions.size() == 1) {
    /* There is a single element in perpendicularPlanePositions. The resulting
     * rotation within that symmetry is the identity rotation, because this
     * single index can be rotated any which way to satisfy the other side.
     */
    return {1};
  }

  /* Remaining case: There are no elements in perpendicularPlanePositions. Then
   * there is no rotation, not even identity, to help in combinatorial handling
   */
  return {};
}

Composite::AngleGroup Composite::smallestAngleGroup(
  const Symmetry::Name symmetryName,
  const unsigned fusedSymmetryPosition
) {
  // Initialize the search state
  AngleGroup angleGroup;
  angleGroup.symmetryPositions.reserve(Symmetry::size(symmetryName));
  angleGroup.angle = M_PI;

  // Go through all symmetry positions excluding the fused symmetry position
  for(unsigned i = 0; i < Symmetry::size(symmetryName); ++i) {
    if(i == fusedSymmetryPosition) {
      continue;
    }

    double angleToFusedPosition = Symmetry::angleFunction(symmetryName)(fusedSymmetryPosition, i);

    // This naturally excludes M_PI angles from the smallest angle group
    if(fpComparator.isLessThan(angleToFusedPosition, angleGroup.angle)) {
      angleGroup.symmetryPositions = {i};
      angleGroup.angle = angleToFusedPosition;
    } else if(fpComparator.isEqual(angleToFusedPosition, angleGroup.angle)) {
      angleGroup.symmetryPositions.push_back(i);
    }
  }

  return angleGroup;
}

Composite::PerpendicularAngleGroups Composite::inGroupAngles(
  const AngleGroup& angleGroup,
  const Symmetry::Name symmetryName
) {
  PerpendicularAngleGroups groups;

  using RecordVector = std::vector<
    std::pair<unsigned, unsigned>
  >;

  temple::forAllPairs(
    angleGroup.symmetryPositions,
    [&](const unsigned a, const unsigned b) -> void {
      double perpendicularAngle = perpendicularSubstituentAngle(
        angleGroup.angle,
        Symmetry::angleFunction(symmetryName)(a, b)
      );

      auto findIter = std::find_if(
        groups.begin(),
        groups.end(),
        [&](const auto& record) -> bool {
          return temple::any_of(
            record.first,
            [&](const double angle) -> bool {
              return fpComparator.isEqual(perpendicularAngle, angle);
            }
          );
        }
      );

      if(findIter == groups.end()) {
        // No equal angles exist, add one yourself
        groups.emplace_back(
          temple::TinyUnorderedSet<double> {perpendicularAngle},
          RecordVector {detail::makeOrderedPair(a, b)}
        );
      } else {
        findIter->second.emplace_back(
          detail::makeOrderedPair(a, b)
        );
      }
    }
  );

  return groups;
}

Composite::Composite(
  const Symmetry::Name left,
  const Symmetry::Name right,
  const unsigned leftFusedPosition,
  const unsigned rightFusedPosition,
  const std::vector<char>& leftCharacters,
  const std::vector<char>& rightCharacters
) : _left(left),
    _right(right),
    _leftFusedPosition(leftFusedPosition),
    _rightFusedPosition(rightFusedPosition)
{
  assert(leftFusedPosition < Symmetry::size(left));
  assert(rightFusedPosition < Symmetry::size(right));

  auto isIsotropic = [](const std::vector<char>& characters) -> bool {
    return temple::all_of(
      characters,
      [&](const char character) -> bool {
        return character == characters.front();
      }
    );
  };

  if(isIsotropic(leftCharacters) || isIsotropic(rightCharacters)) {
    // TODO don't we have to enforce planarity anyway?
    // No permutations possible
    return;
  }

  /* Generate a set of stereopermutations for this particular combination,
   * which can then be indexed
   */

  /* Find the group of symmetry positions with the smallest angle to the
   * fused position (these are the only important ones when considering
   * relative position across both groups).
   */
  auto leftAngleGroup = smallestAngleGroup(left, leftFusedPosition);
  auto rightAngleGroup = smallestAngleGroup(right, rightFusedPosition);

  /*std::cout << "Angle group for " << Symmetry::name(left) << ": "
    << leftAngleGroup.angle << " -> "
    << temple::stringify(leftAngleGroup.symmetryPositions) << "\n";

  std::cout << "Angle group for " << Symmetry::name(right) << ": "
    << rightAngleGroup.angle << " -> "
    << temple::stringify(rightAngleGroup.symmetryPositions) << "\n";*/

  /* Range of possibilities for each side:
   * - no symmetry positions in smallest angle group -> no relative
   *   positioning possible (probably only in linear)
   * - one symmetry position (can be rotated any way it's needed)
   * - multiple symmetry positions (need to determine angle in shared
   *   perpendicular plane to figure out how they can be arranged)
   */

  /* Range of combinatorial possibilities:
   * - Either side has zero symmetry positions in the smallest angle group:
   *   No relative positioning possible, this Composite has zero
   *   stereopermutations
   * - Both sides have one symmetry position in the smallest angle group:
   *   Dihedrals for both can be cis / trans
   * - One side has one symmetry position in the smallest angle group:
   *   Dihedral is 0 to one symmetry position of the larger side, X to the others
   * - Both sides have multiple symmetry positions in the smallest angle group:
   *   Figure out the relative angles between positions in each angle group and
   *   try to find matches across groups -> these can be arranged in a coplanar
   *   fashion. Then each rotation on one side generates a new overlay
   *   possibility.
   */

  /* The central atom of both symmetries is always placed at the origin in the
   * coordinate definitions.
   */
  auto leftCoordinates = Symmetry::symmetryData().at(left).coordinates;
  // Rotate left fused position onto <1, 0, 0>
  detail::rotateCoordinates(
    leftCoordinates,
    leftCoordinates.at(leftFusedPosition).normalized(),
    Eigen::Vector3d::UnitX()
  );

  auto rightCoordinates = Symmetry::symmetryData().at(right).coordinates;
  // Rotate right fused position onto <-1, 0, 0>
  detail::rotateCoordinates(
    rightCoordinates,
    rightCoordinates.at(rightFusedPosition).normalized(),
    -Eigen::Vector3d::UnitX()
  );

  // Translate positions by <1, 0, 0>
  detail::translateCoordinates(
    rightCoordinates,
    Eigen::Vector3d::UnitX()
  );

  auto getDihedral = [&](const unsigned l, const unsigned r) -> double {
    return detail::dihedral(
      leftCoordinates.at(l),
      Eigen::Vector3d::Zero(),
      Eigen::Vector3d::UnitX(),
      rightCoordinates.at(r)
    );
  };

  // TODO abbreviated computation for same symmetry cases?

  /* Sequentially align every pair. Pick that arrangement in which the number
   * of cis dihedrals is maximal.
   *
   * This is essentially brute-forcing the problem. I'm having a hard time
   * thinking up an elegant solution that can satisfy all possible symmetries.
   */
  temple::forAllPairs(
    leftAngleGroup.symmetryPositions,
    rightAngleGroup.symmetryPositions,
    [&](const unsigned l, const unsigned r) -> void {
      // Calculate the dihedral angle from l.front() to r
      double dihedralAngle = getDihedral(l, r);

      // Twist the right coordinates around x so that l.front() is cis with r
      for(auto& position: rightCoordinates) {
        position = Eigen::AngleAxisd(
          -dihedralAngle,
          Eigen::Vector3d::UnitX()
        ) * position;
      }

      // Make sure the rotation leads to cis arrangement
      assert(std::fabs(getDihedral(l, r)) < 1e-10);

      auto dihedralList = temple::mapAllPairs(
        leftAngleGroup.symmetryPositions,
        rightAngleGroup.symmetryPositions,
        [&](const unsigned l, const unsigned r) -> DihedralTuple {
          return {
            l,
            r,
            getDihedral(l, r)
          };
        }
      );

      if(
        !temple::any_of(
          _stereopermutations,
          [&dihedralList](const auto& rhsDihedralList) -> bool {
            return fpComparator.isEqual(
              std::get<2>(dihedralList.front()),
              std::get<2>(rhsDihedralList.front())
            );
          }
        )
      ) {
        _stereopermutations.push_back(std::move(dihedralList));
      }
    }
  );

  /* For situations in which only one position exists in both symmetries, add
   * the trans dihedral possibility explicitly
   */
  if(
    leftAngleGroup.symmetryPositions.size() == 1
    && rightAngleGroup.symmetryPositions.size() == 1
  ) {
    // Add trans dihedral possibility
    _stereopermutations.emplace_back(
      std::vector<DihedralTuple> {
        DihedralTuple {
          leftAngleGroup.symmetryPositions.front(),
          rightAngleGroup.symmetryPositions.front(),
          M_PI
        }
      }
    );
  }

  /*std::cout << "Determined the following " << _stereopermutations.size()
    <<  " stereopermutations:\n";

  for(const auto& permutation : _stereopermutations) {
    std::cout << "- " << temple::stringify(permutation) << "\n";
  }*/
}

unsigned Composite::permutations() const {
  return _stereopermutations.size();
}

const std::vector<Composite::DihedralTuple>& Composite::dihedrals(unsigned permutationIndex) const {
  return _stereopermutations.at(permutationIndex);
}

bool Composite::operator == (const Composite& other) const {
  return (
    _left == other._left
    && _right == other._right
    && _leftFusedPosition == other._leftFusedPosition
    && _rightFusedPosition == other._rightFusedPosition
  );
}

bool Composite::operator != (const Composite& other) const {
  return !(*this == other);
}

} // namespace stereopermutation
