#include "stereopermutation/Composites.h"

#include "boost/dynamic_bitset.hpp"
#include "Eigen/Geometry"

#include "chemical_symmetries/DynamicProperties.h"
#include "chemical_symmetries/Symmetries.h"
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

Composite::OrientationState::OrientationState(
  Symmetry::Name passSymmetry,
  unsigned passFusedPosition,
  std::vector<char> passCharacters,
  std::size_t passIdentifier
) : symmetry(passSymmetry),
    fusedPosition(passFusedPosition),
    characters(std::move(passCharacters)),
    identifier(passIdentifier)
{
  assert(fusedPosition < Symmetry::size(symmetry));
  assert(characters.size() == Symmetry::size(symmetry));
}

void Composite::OrientationState::applyCharacterRotation(
  const std::vector<unsigned>& rotation
) {
  std::vector<char> newCharacters;
  newCharacters.reserve(rotation.size());

  for(const auto& index : rotation) {
    newCharacters.push_back(
      characters.at(index)
    );
  }

  characters = std::move(newCharacters);
}

std::vector<unsigned> Composite::OrientationState::transformToCanonical() {
  /* For canonical comparisons, we must treat all fused positions within the
   * same position group equally. Although the final generated dihedrals must be
   * different (since indexing is still based on the current symmetry positions
   * within each partial symmetry, the sequence must be the same across any
   * position group.
   */
  unsigned reducedFusedPosition = lowestEqualPositionInSymmetry();

  // Find the mapping
  auto toCanonicalMapping = findReductionMapping(reducedFusedPosition);

  // Apply it to the data members of the instance
  fusedPosition = reducedFusedPosition;
  applyCharacterRotation(toCanonicalMapping);

  // Return the inverse mapping to allow reversion to the original state
  return Symmetry::properties::inverseRotation(toCanonicalMapping);
}

void Composite::OrientationState::revert(const std::vector<unsigned>& reversionMapping) {
  // Recover the non-canonical state using the reversion mapping
  applyCharacterRotation(reversionMapping);

  auto findIter = std::find(
    std::begin(reversionMapping),
    std::end(reversionMapping),
    fusedPosition
  );

  assert(findIter != std::end(reversionMapping));

  fusedPosition = findIter - std::begin(reversionMapping);
}

std::vector<unsigned> Composite::OrientationState::findReductionMapping(
  unsigned reducedFusedPosition
) const {
  /* Trivial abbreviation: The identity sequence is viable if the fused
   * position is unchanged. It is the lowest permutation possible, and is hence
   * the solution to this search case.
   */
  if(fusedPosition == reducedFusedPosition) {
    return temple::iota<unsigned>(Symmetry::size(symmetry));
  }

  /* Find a mapping that rotates fusedPosition to reducedFusedPosition. In many
   * cases, there are multiple of these. We can remove these degrees of freedom
   * by choosing that rotation whose resulting permutation has the lowest index
   * of permutation.
   */
  const auto identitySequence = temple::iota<unsigned>(Symmetry::size(symmetry));

  // Track the best rotation
  unsigned lowestIndexOfPermutation = std::numeric_limits<unsigned>::max();
  std::vector<unsigned> bestRotation;

  const unsigned linkLimit = Symmetry::rotations(symmetry).size();
  std::vector<unsigned> chain = {0};
  std::vector<
    std::vector<unsigned>
  > chainRotations = {identitySequence};

  /* Zero-initialize a bitset that can store if an index of permutation has
   * been discovered or not.
   */
  boost::dynamic_bitset<> discoveredIndicesOfPermutation {
    temple::Math::factorial(Symmetry::size(symmetry))
  };

  /* The identity sequence has been discovered (initial element of
   * chainRotations). We set bit 0 since the identity sequence has this index
   * of permutation.
   */
  discoveredIndicesOfPermutation.set(0);

  while(chain.front() < linkLimit) {
    // Generate a new rotation
    auto generatedRotation = Symmetry::properties::applyRotation(
      chainRotations.back(),
      Symmetry::rotations(symmetry).at(
        chain.back()
      )
    );

    unsigned indexOfPermutation = temple::permutationIndex(generatedRotation);
    // Is it new?
    if(discoveredIndicesOfPermutation.test(indexOfPermutation)) {
      // Already discovered! Are we at the maximum instruction?
      if(chain.back() < linkLimit - 1) {
        // No, so we can just increment
        ++chain.back();
      } else {
        /* Collapse the chain until we reach an incrementable position,
         * preserve the first position for the exit condition test
         */
        while(
          chain.size() > 1
          && chain.back() == linkLimit - 1
        ) {
          chain.pop_back();
          chainRotations.pop_back();
        }

        // Increment the back of the chain
        ++chain.back();
      }
    } else {
      // The rotation is new, determine if it should replace the tracked best
      if(
        generatedRotation.at(reducedFusedPosition) == fusedPosition
        && indexOfPermutation < lowestIndexOfPermutation
      ) {
        bestRotation = generatedRotation;
        lowestIndexOfPermutation = indexOfPermutation;
      }

      // Add the new rotation to the set of discovered ones
      discoveredIndicesOfPermutation.set(indexOfPermutation);

      // Add the rotation and a new instruction to the chain
      chainRotations.push_back(std::move(generatedRotation));
      chain.emplace_back(0);
    }
  }

  // We must have found some mapping
  assert(!bestRotation.empty());

  return bestRotation;
}

unsigned Composite::OrientationState::lowestEqualPositionInSymmetry() const {
  auto positionGroupCharacters = Symmetry::properties::positionGroups(symmetry);

  /* Return the position of the first character that matches that of the fused
   * position
   */
  auto findIter = std::find_if(
    std::begin(positionGroupCharacters),
    std::end(positionGroupCharacters),
    [&](const char groupLabel) -> bool {
      return groupLabel == positionGroupCharacters.at(fusedPosition);
    }
  );

  assert(findIter != std::end(positionGroupCharacters));

  return findIter - std::begin(positionGroupCharacters);
}

Composite::AngleGroup Composite::OrientationState::smallestAngleGroup() const {
  // Initialize the search state
  AngleGroup angleGroup;
  angleGroup.symmetryPositions.reserve(Symmetry::size(symmetry));
  angleGroup.angle = M_PI;

  // Go through all symmetry positions excluding the fused symmetry position
  for(unsigned i = 0; i < Symmetry::size(symmetry); ++i) {
    if(i == fusedPosition) {
      continue;
    }

    double angleToFusedPosition = Symmetry::angleFunction(symmetry)(fusedPosition, i);

    // This naturally excludes M_PI angles from the smallest angle group
    if(fpComparator.isLessThan(angleToFusedPosition, angleGroup.angle)) {
      angleGroup.symmetryPositions = {i};
      angleGroup.angle = angleToFusedPosition;
    } else if(fpComparator.isEqual(angleToFusedPosition, angleGroup.angle)) {
      angleGroup.symmetryPositions.push_back(i);
    }
  }

  /* In order to identify if a side is isotropic, check whether the rankings of
   * ligands at this angle group are all the same.
   */
  if(angleGroup.symmetryPositions.size() == 1) {
    // A single relevant symmetry position is not isotropic
    angleGroup.isotropic = false;
  } else {
    auto relevantCharacters = temple::map(
      angleGroup.symmetryPositions,
      [&](const unsigned symmetryPosition) -> char {
        return characters.at(symmetryPosition);
      }
    );

    angleGroup.isotropic = temple::all_of(
      relevantCharacters,
      [&](const char character) -> bool {
        return character == relevantCharacters.front();
      }
    );
  }

  return angleGroup;
}

bool Composite::OrientationState::operator < (const Composite::OrientationState& other) const {
  return (
    std::tie(symmetry, fusedPosition, characters)
    < std::tie(other.symmetry, other.fusedPosition, other.characters)
  );
}

bool Composite::OrientationState::operator == (const Composite::OrientationState& other) const {
  return (
    std::tie(symmetry, fusedPosition, characters)
    == std::tie(other.symmetry, other.fusedPosition, other.characters)
  );
}

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
      rotation = temple::iota<unsigned>(Symmetry::size(symmetryName));

      for(const auto r : rotationIndexApplicationSequence) {
        rotation = Symmetry::properties::applyRotation(
          rotation,
          symmetryName,
          r
        );
      }

      // Determine if the current rotation matches all required criteria
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

Composite::Composite(OrientationState first, OrientationState second)
  : _orientations { std::move(first), std::move(second) }
{
  // Do not construct the ordered pair of OrientationStates with same identifier
  assert(_orientations.first.identifier != _orientations.second.identifier);

  /* For canonical comparison, the OrientationStates must be transformed as
   * though the fused position was the lowest position within the symmetry
   * position group, and then back-transformed later.
   */

  auto firstReversionMapping = _orientations.first.transformToCanonical();
  auto secondReversionMapping = _orientations.second.transformToCanonical();

  /* Find the group of symmetry positions with the smallest angle to the
   * fused position (these are the only important ones when considering
   * relative position across both groups).
   */
  auto angleGroups = _orientations.map(
    [](const OrientationState& orientation) -> AngleGroup {
      return orientation.smallestAngleGroup();
    }
  );

  /* Even if either side is isotropic and overall this stereocenter has only one
   * assignment, the first discovered permutation must be kept in order to
   * enforce planarity if it exists. Any further stereopermutations lead to
   * identical three-dimensional arrangements.
   *
   */
  bool overallIsotropic = angleGroups.first.isotropic || angleGroups.second.isotropic;

  /* Generate a set of stereopermutations for this particular combination,
   * which can then be indexed
   *
   * Range of combinatorial possibilities:
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
   *
   * NOTE: The central atom of both symmetries is always placed at the origin
   * in the coordinate definitions.
   */
  auto firstCoordinates = Symmetry::symmetryData().at(
    _orientations.first.symmetry
  ).coordinates;
  // Rotate left fused position onto <1, 0, 0>
  detail::rotateCoordinates(
    firstCoordinates,
    firstCoordinates.at(_orientations.first.fusedPosition).normalized(),
    Eigen::Vector3d::UnitX()
  );

  auto secondCoordinates = Symmetry::symmetryData().at(
    _orientations.second.symmetry
  ).coordinates;
  // Rotate right fused position onto <-1, 0, 0>
  detail::rotateCoordinates(
    secondCoordinates,
    secondCoordinates.at(_orientations.second.fusedPosition).normalized(),
    -Eigen::Vector3d::UnitX()
  );

  // Translate positions by <1, 0, 0>
  detail::translateCoordinates(
    secondCoordinates,
    Eigen::Vector3d::UnitX()
  );

  auto getDihedral = [&](const unsigned f, const unsigned s) -> double {
    return detail::dihedral(
      firstCoordinates.at(f),
      Eigen::Vector3d::Zero(),
      Eigen::Vector3d::UnitX(),
      secondCoordinates.at(s)
    );
  };

  // TODO is there a simple abbreviated computation for same symmetry cases?

  /* Sequentially align every pair. Pick that arrangement in which the number
   * of cis dihedrals is maximal.
   *
   * This is essentially brute-forcing the problem. I'm having a hard time
   * thinking up an elegant solution that can satisfy all possible symmetries.
   */
  auto addDihedralCombination = [&](const unsigned f, const unsigned s) -> void {
    // Calculate the dihedral angle from l.front() to r
    double dihedralAngle = getDihedral(f, s);

    // Twist the right coordinates around x so that l.front() is cis with r
    for(auto& position: secondCoordinates) {
      position = Eigen::AngleAxisd(
        -dihedralAngle,
        Eigen::Vector3d::UnitX()
      ) * position;
    }

    // Make sure the rotation leads to cis arrangement
    assert(std::fabs(getDihedral(f, s)) < 1e-10);

    auto dihedralList = temple::mapAllPairs(
      angleGroups.first.symmetryPositions,
      angleGroups.second.symmetryPositions,
      [&](const unsigned f, const unsigned s) -> DihedralTuple {
        return {
          f,
          s,
          getDihedral(f, s)
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
  };

  // Branch if only one stereopermutation is needed due to isotropicity
  if(overallIsotropic) {
    addDihedralCombination(
      angleGroups.first.symmetryPositions.front(),
      angleGroups.second.symmetryPositions.front()
    );
  } else {
    temple::forAllPairs(
      angleGroups.first.symmetryPositions,
      angleGroups.second.symmetryPositions,
      addDihedralCombination
    );
  }

  /* For situations in which only one position exists in both symmetries, add
   * the trans dihedral possibility explicitly
   */
  if(
    angleGroups.first.symmetryPositions.size() == 1
    && angleGroups.second.symmetryPositions.size() == 1
  ) {
    // Add trans dihedral possibility
    _stereopermutations.emplace_back(
      std::vector<DihedralTuple> {
        DihedralTuple {
          angleGroups.first.symmetryPositions.front(),
          angleGroups.second.symmetryPositions.front(),
          M_PI
        }
      }
    );
  }

  // Revert the OrientationStates and transform the stereopermutations too
  _orientations.first.revert(firstReversionMapping);
  _orientations.second.revert(secondReversionMapping);

  for(auto& stereopermutation : _stereopermutations) {
    for(auto& dihedralTuple : stereopermutation) {
      std::get<0>(dihedralTuple) = firstReversionMapping.at(
        std::get<0>(dihedralTuple)
      );

      std::get<1>(dihedralTuple) = secondReversionMapping.at(
        std::get<1>(dihedralTuple)
      );
    }
  }

}

unsigned Composite::permutations() const {
  return _stereopermutations.size();
}

const std::vector<Composite::DihedralTuple>& Composite::dihedrals(unsigned permutationIndex) const {
  return _stereopermutations.at(permutationIndex);
}

const temple::OrderedPair<Composite::OrientationState>& Composite::orientations() const {
  return _orientations;
}

Composite::PermutationsList::const_iterator Composite::begin() const {
  return std::begin(_stereopermutations);
}

Composite::PermutationsList::const_iterator Composite::end() const {
  return std::end(_stereopermutations);
}

bool Composite::operator == (const Composite& other) const {
  return _orientations == other._orientations;
}

bool Composite::operator != (const Composite& other) const {
  return !(*this == other);
}

} // namespace stereopermutation
