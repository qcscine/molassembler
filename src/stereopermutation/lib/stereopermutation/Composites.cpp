/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "stereopermutation/Composites.h"

#include "boost/dynamic_bitset.hpp"
#include "Eigen/Geometry"

#include "shapes/Properties.h"
#include "shapes/Data.h"
#include "temple/Adaptors/AllPairs.h"
#include "temple/Functor.h"
#include "temple/Functional.h"
#include "temple/Permutations.h"
#include "temple/Stringify.h"

namespace Scine {

namespace stereopermutation {

namespace detail {

template<typename T>
std::pair<T, T> makeOrderedPair(T a, T b) {
  std::pair<T, T> pair {
    std::move(a),
    std::move(b)
  };

  if(pair.second < pair.first) {
    std::swap(pair.first, pair.second);
  }

  return pair;
}

void rotateCoordinates(
  Eigen::Ref<shapes::CoordinateList> positions,
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
    positions *= -1;
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

  for(unsigned i = 0; i < positions.cols(); ++i) {
    positions.col(i) = rotation * positions.col(i);
  }
}

void translateCoordinates(
  Eigen::Ref<shapes::CoordinateList> positions,
  const Eigen::Vector3d& translation
) {
  for(unsigned i = 0; i < positions.cols(); ++i) {
    positions.col(i) += translation;
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
  shapes::Shape passShape,
  unsigned passFusedVertex,
  std::vector<char> passCharacters,
  std::size_t passIdentifier
) : shape(passShape),
    fusedVertex(passFusedVertex),
    characters(std::move(passCharacters)),
    identifier(passIdentifier)
{
  assert(fusedVertex < shapes::size(shape));
  assert(characters.size() == shapes::size(shape));
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
   * different (since indexing is still based on the current shape positions
   * within each partial shape, the sequence must be the same across any
   * position group.
   */
  unsigned reducedFusedVertex = lowestEqualVertexInShape();

  // Find the mapping
  auto toCanonicalMapping = findReductionMapping(reducedFusedVertex);

  // Apply it to the data members of the instance
  fusedVertex = reducedFusedVertex;
  applyCharacterRotation(toCanonicalMapping);

  // Return the inverse mapping to allow reversion to the original state
  return shapes::properties::inverseRotation(toCanonicalMapping);
}

void Composite::OrientationState::revert(const std::vector<unsigned>& reversionMapping) {
  // Recover the non-canonical state using the reversion mapping
  applyCharacterRotation(reversionMapping);

  auto findIter = std::find(
    std::begin(reversionMapping),
    std::end(reversionMapping),
    fusedVertex
  );

  assert(findIter != std::end(reversionMapping));

  fusedVertex = findIter - std::begin(reversionMapping);
}

std::vector<unsigned> Composite::OrientationState::findReductionMapping(
  unsigned reducedFusedVertex
) const {
  /* NOTE: The implementation below is VERY similar to
   * Stereopermutation::_generateAllRotation's generation work.
   * BUT! This algorithm doesn't store the rotations, merely their index of
   * permutation to be able to terminate the backtracking algorithm, and tracks
   * the best structure.
   *
   * Refactoring both to a common denominator could be challenging.
   */

  /* Trivial abbreviation: The identity sequence is viable if the fused
   * position is unchanged. It is the lowest permutation possible, and is hence
   * the solution to this search case.
   */
  if(fusedVertex == reducedFusedVertex) {
    return temple::iota<unsigned>(shapes::size(shape));
  }

  /* Find a mapping that rotates fusedVertex to reducedFusedVertex. In many
   * cases, there are multiple of these. We can remove these degrees of freedom
   * by choosing that rotation whose resulting permutation has the lowest index
   * of permutation.
   */
  const auto identitySequence = temple::iota<unsigned>(shapes::size(shape));

  // Track the best rotation
  unsigned lowestIndexOfPermutation = std::numeric_limits<unsigned>::max();
  std::vector<unsigned> bestRotation;

  const unsigned linkLimit = shapes::rotations(shape).size();
  std::vector<unsigned> chain = {0};
  chain.reserve(32);
  std::vector<
    std::vector<unsigned>
  > chainRotations = {identitySequence};

  /* Zero-initialize a bitset that can store if an index of permutation has
   * been discovered or not.
   */
  boost::dynamic_bitset<> discoveredIndicesOfPermutation {
    temple::Math::factorial(shapes::size(shape))
  };

  /* The identity sequence has been discovered (initial element of
   * chainRotations). We set bit 0 since the identity sequence has this index
   * of permutation.
   */
  discoveredIndicesOfPermutation.set(0);

  while(chain.front() < linkLimit) {
    // Generate a new rotation
    auto generatedRotation = shapes::properties::applyRotation(
      chainRotations.back(),
      shapes::rotations(shape).at(
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
        generatedRotation.at(reducedFusedVertex) == fusedVertex
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

unsigned Composite::OrientationState::lowestEqualVertexInShape() const {
  const auto positionGroupCharacters = shapes::properties::positionGroups(shape);

  /* Return the position of the first character that matches that of the fused
   * position
   */
  auto findIter = std::find_if(
    std::begin(positionGroupCharacters),
    std::end(positionGroupCharacters),
    [&](const char groupLabel) -> bool {
      return groupLabel == positionGroupCharacters.at(fusedVertex);
    }
  );

  assert(findIter != std::end(positionGroupCharacters));

  return findIter - std::begin(positionGroupCharacters);
}

Composite::AngleGroup Composite::OrientationState::smallestAngleGroup() const {
  // Initialize the search state
  AngleGroup angleGroup;
  angleGroup.shapeVertices.reserve(shapes::size(shape));
  angleGroup.angle = M_PI;

  // Go through all symmetry positions excluding the fused shape position
  for(unsigned i = 0; i < shapes::size(shape); ++i) {
    if(i == fusedVertex) {
      continue;
    }

    double angleToFusedPosition = shapes::angleFunction(shape)(fusedVertex, i);

    // This naturally excludes M_PI angles from the smallest angle group
    if(fpComparator.isLessThan(angleToFusedPosition, angleGroup.angle)) {
      angleGroup.shapeVertices = {i};
      angleGroup.angle = angleToFusedPosition;
    } else if(fpComparator.isEqual(angleToFusedPosition, angleGroup.angle)) {
      angleGroup.shapeVertices.push_back(i);
    }
  }

  /* In order to identify if a side is isotropic, check whether the rankings of
   * ligands at this angle group are all the same.
   */
  if(angleGroup.shapeVertices.size() == 1) {
    // A single relevant symmetry position is not isotropic
    angleGroup.isotropic = false;
  } else {
    auto relevantCharacters = temple::map(
      angleGroup.shapeVertices,
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

  assert(
    std::is_sorted(
      std::begin(angleGroup.shapeVertices),
      std::end(angleGroup.shapeVertices)
    )
  );

  return angleGroup;
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
  const shapes::Shape shape,
  const unsigned fixedSymmetryPosition,
  const std::vector<unsigned>& changedPositions
) {
  auto periodicities = temple::map(
    temple::iota<unsigned>(shapes::rotations(shape).size()),
    [&shape](const unsigned rotationFunctionIndex) -> unsigned {
      return shapes::properties::rotationPeriodicity(
        shape,
        shapes::rotations(shape).at(rotationFunctionIndex)
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

    temple::inplace::sort(rotationIndexApplicationSequence);

    do {
      // Create the rotation using the index application sequence front-to-back
      rotation = temple::iota<unsigned>(shapes::size(shape));

      for(const auto r : rotationIndexApplicationSequence) {
        rotation = shapes::properties::applyRotation(
          rotation,
          shape,
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
      && temple::inplace::next_permutation(rotationIndexApplicationSequence)
    );

  } while(
    !rotationFound
    && temple::inplace::nextCombinationPermutation(rotationUses, periodicities)
  );

  if(rotationFound) {
    return rotation;
  }

  return {};
}

std::vector<unsigned> Composite::rotation(
  const shapes::Shape shape,
  const unsigned fixedSymmetryPosition,
  const std::vector<unsigned>& perpendicularPlanePositions
) {
  // Three possibilities:

  if(perpendicularPlanePositions.size() > 1) {
    /* There are multiple elements in perpendicularPlanePositions. We have to
     * generate a rotation that keeps fixedSymmetryPosition fixed but rotates the
     * perpendicularPlanePositions, ideally with a periodicity equivalent to the
     * amount of symmetry positions involved.
     */
    auto candidateRotation = generateRotation(
      shape,
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
      shapes::properties::rotationPeriodicity(shape, candidateRotation)
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
  const shapes::Shape shape
) {
  PerpendicularAngleGroups groups;

  using RecordVector = std::vector<
    std::pair<unsigned, unsigned>
  >;

  temple::forEach(
    temple::adaptors::allPairs(angleGroup.shapeVertices),
    [&](const unsigned a, const unsigned b) -> void {
      double perpendicularAngle = perpendicularSubstituentAngle(
        angleGroup.angle,
        shapes::angleFunction(shape)(a, b)
      );

      auto findIter = std::find_if(
        std::begin(groups),
        std::end(groups),
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
          std::vector<double> {perpendicularAngle},
          RecordVector {detail::makeOrderedPair(a, b)}
        );
      } else {
        findIter->first.push_back(perpendicularAngle);
        findIter->second.emplace_back(
          detail::makeOrderedPair(a, b)
        );
      }
    }
  );

  return groups;
}

Composite::Composite(
  OrientationState first,
  OrientationState second,
  const Alignment alignment
) : _orientations { std::move(first), std::move(second) },
    _alignment(alignment)
{
  // Do not construct the ordered pair of OrientationStates with same identifier
  assert(_orientations.first.identifier != _orientations.second.identifier);

  /* In order to get meaningful indices of permutation, combinations of
   * symmetries across fused positions within the same group of symmetry
   * positions (e.g. equatorial or apical in square pyramidal) must be the same.
   *
   * In order to achieve this, the OrientationStates is transformed by a
   * rotation that temporarily places the fused position at the lowest index
   * symmetry position in its symmetry. After permutation are generated,
   * the orientation state is transformed back.
   */
  auto firstReversionMapping = _orientations.first.transformToCanonical();
  auto secondReversionMapping = _orientations.second.transformToCanonical();

  /* Find the group of symmetry positions with the smallest angle to the
   * fused position (these are the only important ones when considering
   * relative arrangements across the bond).
   */
  auto angleGroups = _orientations.map(
    [](const OrientationState& orientation) -> AngleGroup {
      return orientation.smallestAngleGroup();
    }
  );

  /* Reorder both AngleGroups' shapeVertices by descending ranking and
   * index to get canonical initial combinations
   */
  temple::inplace::sort(
    angleGroups.first.shapeVertices,
    [&](const unsigned a, const unsigned b) -> bool {
      return (
        std::tie(_orientations.first.characters.at(a), a)
        > std::tie(_orientations.first.characters.at(b), b)
      );
    }
  );

  temple::inplace::sort(
    angleGroups.second.shapeVertices,
    [&](const unsigned a, const unsigned b) -> bool {
      return (
        std::tie(_orientations.second.characters.at(a), a)
        > std::tie(_orientations.second.characters.at(b), b)
      );
    }
  );

  /* From the angle groups' characters, we can figure out if all
   * stereopermutations that are generated will be ranking-wise equivalent
   * spatially despite differing in the symmetry positions at which the equally
   * ranked substituents are placed. This is important information for deciding
   * whether a Composite yields a stereogenic object.
   */
  _isotropic = angleGroups.first.isotropic || angleGroups.second.isotropic;

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
  shapes::CoordinateList firstCoordinates = shapes::shapeData().at(
    _orientations.first.shape
  ).coordinates;
  // Rotate left fused position onto <1, 0, 0>
  detail::rotateCoordinates(
    firstCoordinates,
    firstCoordinates.col(_orientations.first.fusedVertex).normalized(),
    Eigen::Vector3d::UnitX()
  );

  auto secondCoordinates = shapes::shapeData().at(
    _orientations.second.shape
  ).coordinates;
  // Rotate right fused position onto <-1, 0, 0>
  detail::rotateCoordinates(
    secondCoordinates,
    secondCoordinates.col(_orientations.second.fusedVertex).normalized(),
    -Eigen::Vector3d::UnitX()
  );

  // Translate positions by <1, 0, 0>
  detail::translateCoordinates(
    secondCoordinates,
    Eigen::Vector3d::UnitX()
  );

  auto getDihedral = [&](const unsigned f, const unsigned s) -> double {
    return detail::dihedral(
      firstCoordinates.col(f),
      Eigen::Vector3d::Zero(),
      Eigen::Vector3d::UnitX(),
      secondCoordinates.col(s)
    );
  };

  /* Sequentially align every pair. Pick that arrangement in which the number
   * of cis dihedrals is maximal.
   *
   * This is essentially brute-forcing the problem. I'm having a hard time
   * thinking up an elegant solution that can satisfy all possible symmetries.
   */

  // Generate all arrangements regardless of whether the Composite is isotropic
  temple::forEach(
    temple::adaptors::allPairs(
      angleGroups.first.shapeVertices,
      angleGroups.second.shapeVertices
    ),
    [&](const unsigned f, const unsigned s) -> void {
      const double alignAngle = getDihedral(f, s);

      // Twist the right coordinates around x so that f is cis with r
      for(unsigned i = 0; i < secondCoordinates.cols(); ++i) {
        secondCoordinates.col(i) = Eigen::AngleAxisd(
          -alignAngle,
          Eigen::Vector3d::UnitX()
        ) * secondCoordinates.col(i);
      }

      // Make sure the rotation leads to a zero dihedral
      assert(std::fabs(getDihedral(f, s)) < 1e-10);

      // Offset if desired
      double offsetAngle = 0;
      if(alignment == Alignment::Staggered) {
        /* The offset angle for a staggered arrangement is half of the angle
         * to the next symmetry position in some direction. We'll go for most
         * negative.
         */

        const auto dihedrals = temple::map(
          angleGroups.second.shapeVertices,
          [&](const unsigned secondSymmetryPosition) -> double {
            double dihedral = getDihedral(f, secondSymmetryPosition);
            if(dihedral >= -1e-10) {
              dihedral -= 2 * M_PI;
            }
            return dihedral;
          }
        );

        const unsigned maximumIndex = std::max_element(
          std::begin(dihedrals),
          std::end(dihedrals)
        ) - std::begin(dihedrals);

        // std::cout << "Next symmetry position in rotor is " << angleGroups.second.shapeVertices.at(maximumIndex) << " with dihedral of " << dihedrals.at(maximumIndex) << "\n";

        offsetAngle = dihedrals.at(maximumIndex) / 2;
      }

      if(offsetAngle != 0.0) {
        for(unsigned i = 0; i < secondCoordinates.cols(); ++i) {
          secondCoordinates.col(i) = Eigen::AngleAxisd(
            offsetAngle,
            Eigen::Vector3d::UnitX()
          ) * secondCoordinates.col(i);
        }
      }

      auto dihedralList = temple::map(
        temple::adaptors::allPairs(
          angleGroups.first.shapeVertices,
          angleGroups.second.shapeVertices
        ),
        [&](const unsigned a, const unsigned b) -> DihedralTuple {
          return {a, b, getDihedral(a, b)};
        }
      );

      // Ensure postcondition that list of dihedrals is sorted
      std::sort(
        std::begin(dihedralList),
        std::end(dihedralList)
      );

      if(
        !temple::any_of(
          _stereopermutations,
          [&dihedralList](const auto& rhsDihedralList) -> bool {
            return fpComparator.isEqual(
              std::get<2>(dihedralList.front()),
              std::get<2>(rhsDihedralList.front())
            ) || (
              fpComparator.isEqual(
                std::fabs(std::get<2>(dihedralList.front())),
                M_PI
              ) && fpComparator.isEqual(
                std::fabs(std::get<2>(rhsDihedralList.front())),
                M_PI
              )
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
    angleGroups.first.shapeVertices.size() == 1
    && angleGroups.second.shapeVertices.size() == 1
  ) {
    // Add trans dihedral possibility
    _stereopermutations.emplace_back(
      std::vector<DihedralTuple> {
        DihedralTuple {
          angleGroups.first.shapeVertices.front(),
          angleGroups.second.shapeVertices.front(),
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
      auto findIter = std::find(
        std::begin(firstReversionMapping),
        std::end(firstReversionMapping),
        std::get<0>(dihedralTuple)
      );

      assert(findIter != std::end(firstReversionMapping));

      std::get<0>(dihedralTuple) = findIter - std::begin(firstReversionMapping);

      assert(std::get<0>(dihedralTuple) != _orientations.first.fusedVertex);

      findIter = std::find(
        std::begin(secondReversionMapping),
        std::end(secondReversionMapping),
        std::get<1>(dihedralTuple)
      );

      assert(findIter != std::end(secondReversionMapping));

      std::get<1>(dihedralTuple) = findIter - std::begin(secondReversionMapping);

      assert(std::get<1>(dihedralTuple) != _orientations.second.fusedVertex);
    }
  }

  /* Reverse the stereopermutation sequence. This is so that the indices of the
   * generated permutations yield the following simple comparison:
   *
   *   0 is E, 1 is Z
   *   1 > 0 == Z > E
   */
  std::reverse(
    std::begin(_stereopermutations),
    std::end(_stereopermutations)
  );
}

void Composite::applyIdentifierPermutation(const std::vector<std::size_t>& permutation) {
  for(auto& orientationState : _orientations) {
    orientationState.identifier = permutation.at(orientationState.identifier);
  }
}

unsigned Composite::permutations() const {
  return _stereopermutations.size();
}

Composite::Alignment Composite::alignment() const {
  return _alignment;
}

const std::vector<Composite::DihedralTuple>& Composite::dihedrals(unsigned permutationIndex) const {
  return _stereopermutations.at(permutationIndex);
}

bool Composite::isIsotropic() const {
  return _isotropic;
}

unsigned Composite::order() const {
  auto countDistinct = [&](auto&& f) {
    std::set<unsigned> positions;
    for(const DihedralTuple& t : _stereopermutations.front()) {
      positions.insert(f(t));
    }
    return positions.size();
  };

  return std::max(
    countDistinct(temple::functor::first),
    countDistinct(temple::functor::second)
  );
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

bool Composite::operator < (const Composite& other) const {
  return _orientations < other._orientations;
}

bool Composite::operator == (const Composite& other) const {
  return _orientations == other._orientations;
}

bool Composite::operator != (const Composite& other) const {
  return !(*this == other);
}

} // namespace stereopermutation

} // namespace Scine
