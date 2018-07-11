#include "Composites.h"

#include "chemical_symmetries/DynamicProperties.h"
#include "temple/Containers.h"
#include "temple/Stringify.h"

#include <iostream>

namespace stereopermutation {

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
  const std::vector<unsigned> changedPositions
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

  auto rotationAltersPositions = [&](const std::vector<unsigned> rotation) -> bool {
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
      return identityRotation(symmetryName);
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
    return identityRotation(symmetryName);
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
  const unsigned rightFusedPosition
) : _left(left),
    _right(right),
    _leftFusedPosition(leftFusedPosition),
    _rightFusedPosition(rightFusedPosition)
{
  assert(leftFusedPosition < Symmetry::size(left));
  assert(rightFusedPosition < Symmetry::size(right));

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

  /* Instead of trying to gain understanding of the complete picture, perhaps
   * a minimal solution is sufficient where dihedrals for one index of one side
   * to all others is specified. The remaining information should be specified
   * completely by the two stereocenters on either end of the composite and / or
   * be inferable with triangle inequalities.
   */

  if(leftAngleGroup.symmetryPositions.size() <= rightAngleGroup.symmetryPositions.size()) {
    // Handle the case for left has fewer
    unsigned pickedIndex = leftAngleGroup.symmetryPositions.front();

    _stereopermutations.resize(rightAngleGroup.symmetryPositions.size());

    for(unsigned i = 0; i < rightAngleGroup.symmetryPositions.size(); ++i) {
      auto& dihedrals = _stereopermutations.at(i);

      // Model the cis dihedral to right's symmetryPosition i
      dihedrals.emplace_back(
        pickedIndex,
        rightAngleGroup.symmetryPositions.at(i),
        0.0
      );

      // And to all others in right if picked - i is cis
      for(unsigned j = 0; j < rightAngleGroup.symmetryPositions.size(); ++j) {
        if(i == j) { // Skip identical indices case
          continue;
        }

        dihedrals.emplace_back(
          pickedIndex,
          rightAngleGroup.symmetryPositions.at(j),
          perpendicularSubstituentAngle(
            rightAngleGroup.angle,
            Symmetry::angleFunction(right)(
              rightAngleGroup.symmetryPositions.at(i),
              rightAngleGroup.symmetryPositions.at(j)
            )
          )
        );
      }
    }
  } else {
    // In case right has fewer
    unsigned pickedIndex = rightAngleGroup.symmetryPositions.front();

    _stereopermutations.resize(leftAngleGroup.symmetryPositions.size());

    for(unsigned i = 0; i < leftAngleGroup.symmetryPositions.size(); ++i) {
      auto& dihedrals = _stereopermutations.at(i);

      // Model the cis dihedral to left's symmetryPosition i
      dihedrals.emplace_back(
        leftAngleGroup.symmetryPositions.at(i),
        pickedIndex,
        0.0
      );

      // And to all others in left if picked - i is cis
      for(unsigned j = 0; j < leftAngleGroup.symmetryPositions.size(); ++j) {
        if(i == j) { // Skip identical indices case
          continue;
        }

        dihedrals.emplace_back(
          leftAngleGroup.symmetryPositions.at(j),
          pickedIndex,
          perpendicularSubstituentAngle(
            leftAngleGroup.angle,
            Symmetry::angleFunction(left)(
              leftAngleGroup.symmetryPositions.at(i),
              leftAngleGroup.symmetryPositions.at(j)
            )
          )
        );
      }
    }
  }

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

  /*std::cout << "Determined the following stereopermutations : "
    << temple::stringify(_stereopermutations) << "\n";*/
}

unsigned Composite::permutations() const {
  return _stereopermutations.size();
}

const std::vector<Composite::DihedralTuple>& Composite::dihedrals(unsigned permutationIndex) const {
  return _stereopermutations.at(permutationIndex);
}

/* --------------------------------------------------------------------
 * --------------------------------------------------------------------
 * Dead code for reference if needed later for a more complete solution
 * --------------------------------------------------------------------
 * --------------------------------------------------------------------
 */

//  /* Look for rotations that keep the fused symmetry position constant but
//   * move the symmetry positions in the respective angle groups.
//   */
//  auto leftRotation = rotation(left, leftFusedPosition, leftAngleGroup.symmetryPositions);
//  auto rightRotation = rotation(right, rightFusedPosition, rightAngleGroup.symmetryPositions);
//
//  std::cout << "Left rotation: " << temple::stringify(leftRotation) << "\n";
//  std::cout << "Right rotation: " << temple::stringify(rightRotation) << "\n";
//
//  auto leftInGroupAngles = inGroupAngles(leftAngleGroup, left);
//  auto rightInGroupAngles = inGroupAngles(rightAngleGroup, right);
//
//  std::cout << "In group angle pairs for " << Symmetry::name(left) << ":\n";
//  for(const auto& record : leftInGroupAngles) {
//    std::cout << temple::stringify(record.first.data) << " -> "
//      << temple::stringify(record.second) << "\n";
//  }
//
//  std::cout << "In group angle pairs for " << Symmetry::name(right) << ":\n";
//  for(const auto& record : rightInGroupAngles) {
//    std::cout << temple::stringify(record.first.data) << " -> "
//      << temple::stringify(record.second) << "\n";
//  }
//
//  // Now look for matching angle sets across both symmetries
//  std::vector<
//    std::pair<
//      std::vector<
//        std::pair<unsigned, unsigned>
//      >,
//      std::vector<
//        std::pair<unsigned, unsigned>
//      >
//    >
//  > matches;
//
//  temple::forAllPairs(
//    leftInGroupAngles,
//    rightInGroupAngles,
//    [&](const auto& leftRecord, const auto& rightRecord) -> void {
//      if(
//        temple::any_of(
//          temple::mapAllPairs(
//            leftRecord.first,
//            rightRecord.first,
//            [&](const double leftAngle, const double rightAngle) -> bool {
//              return fpComparator.isEqual(leftAngle, rightAngle);
//            }
//          )
//        )
//      ) {
//        matches.push_back(
//          std::make_pair(leftRecord.second, rightRecord.second)
//        );
//      }
//    }
//  );
//
//  /* There can be multiple matches, e.g. in octahedral - octahedral, where
//   * there are sets of 90° and 180° in each. For those ambiguous cases, we
//   * select the match where the most symmetry positions are involved.
//   */
//  auto matchIter = std::max_element(
//    std::begin(matches),
//    std::end(matches),
//    [](const auto& a, const auto& b) -> unsigned {
//      return (
//        a.first.size() + a.second.size()
//        < b.first.size() + b.second.size()
//      );
//    }
//  );
//
//  /* matchIter is matches.end() if
//   * - Either symmetry has only one symmetry position in angleGroup
//   */
//
//  /*if(matchIter == matches.end()) {
//    std::cout << "No match between symmetry positions of "
//      << Symmetry::name(left) << " and " << Symmetry::name(right) << " found.\n";
//  } else {
//    std::cout << "Selected match for "
//      << Symmetry::name(left) << " and " << Symmetry::name(right) << ": "
//      << temple::stringify(*matchIter) << "\n";
//  }*/

} // namespace stereopermutation
