/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutation/Composites.h"

#include "boost/dynamic_bitset.hpp"
#include "Eigen/Geometry"

#include "Molassembler/Shapes/Properties.h"
#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Functor.h"
#include "Molassembler/Temple/GroupBy.h"
#include "Molassembler/Temple/Permutations.h"

namespace Scine {
namespace Molassembler {
namespace Stereopermutations {
namespace Detail {

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
  Eigen::Ref<Shapes::Coordinates> positions,
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
  const Eigen::Vector3d v = unitSource.cross(unitTarget);

  // Skew-symmetric cross product matrix
  Eigen::Matrix3d v_x;
  v_x <<     0, -v.z(),  v.y(),
         v.z(),      0, -v.x(),
        -v.y(),  v.x(),      0;

  // Dot product (essentially the cosine of the angle for these unit vectors)
  const double c = unitSource.dot(unitTarget);

  // Calculate the rotation matrix
  Eigen::Matrix3d rotation;
  rotation = Eigen::Matrix3d::Identity() + v_x + v_x * v_x * (1.0 / (1 + c));

  for(unsigned i = 0; i < positions.cols(); ++i) {
    positions.col(i) = rotation * positions.col(i);
  }
}

void translateCoordinates(
  Eigen::Ref<Shapes::Coordinates> positions,
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
  const Eigen::Vector3d a = j - i;
  const Eigen::Vector3d b = k - j;
  const Eigen::Vector3d c = l - k;

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

bool dihedralClose(
  const double a,
  const double b,
  const double epsilon
) {
  double diff = b - a;

  if(diff > M_PI) {
    diff -= 2 * M_PI;
  }

  if(diff <= -M_PI) {
    diff += 2 * M_PI;
  }

  return std::fabs(diff) < epsilon;
}

} // namespace Detail

constexpr Temple::Floating::ExpandedAbsoluteEqualityComparator<double> Composite::fpComparator;

Composite::OrientationState::OrientationState(
  Shapes::Shape passShape,
  Shapes::Vertex passFusedVertex,
  std::vector<char> passCharacters,
  std::size_t passIdentifier
) : shape(passShape),
    fusedVertex(passFusedVertex),
    characters(std::move(passCharacters)),
    identifier(passIdentifier)
{
  assert(fusedVertex < Shapes::size(shape));
  assert(characters.size() == Shapes::size(shape));
}

std::vector<char> Composite::OrientationState::applyCharacterRotation(
  const std::vector<Shapes::Vertex>& rotation
) const {
  std::vector<char> newCharacters;
  newCharacters.reserve(rotation.size());

  for(const auto& index : rotation) {
    newCharacters.push_back(
      characters.at(index)
    );
  }

  return newCharacters;
}

std::vector<Shapes::Vertex> Composite::OrientationState::transformToCanonical() {
  /* For canonical comparisons, we must treat all fused positions within the
   * same position group equally. Although the final generated dihedrals must be
   * different (since indexing is still based on the current shape positions
   * within each partial shape, the sequence must be the same across any
   * position group.
   */
  Shapes::Vertex reducedFusedVertex = lowestEqualVertexInShape();

  // Find the mapping
  auto toCanonicalMapping = findReductionMapping(reducedFusedVertex);

  // Apply it to the data members of the instance
  fusedVertex = reducedFusedVertex;
  characters = applyCharacterRotation(toCanonicalMapping);

  // Return the inverse mapping to allow reversion to the original state
  return Shapes::Properties::inverseRotation(toCanonicalMapping);
}

void Composite::OrientationState::revert(const std::vector<Shapes::Vertex>& reversionMapping) {
  // Recover the non-canonical state using the reversion mapping
  characters = applyCharacterRotation(reversionMapping);

  auto findIter = std::find(
    std::begin(reversionMapping),
    std::end(reversionMapping),
    fusedVertex
  );

  assert(findIter != std::end(reversionMapping));

  fusedVertex = findIter - std::begin(reversionMapping);
}

std::vector<Shapes::Vertex> Composite::OrientationState::findReductionMapping(
  const Shapes::Vertex reducedFusedVertex
) const {
  /* NOTE: The implementation below is VERY similar to
   * Stereopermutation::generateAllRotation_'s generation work.
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
    return Temple::iota<Shapes::Vertex>(Shapes::size(shape));
  }

  /* Find a mapping that rotates fusedVertex to reducedFusedVertex. In many
   * cases, there are multiple of these. We can remove these degrees of freedom
   * by choosing that rotation whose resulting permutation has the lowest index
   * of permutation.
   */
  const auto identitySequence = Temple::iota<Shapes::Vertex>(Shapes::size(shape));

  // Track the best rotation
  unsigned lowestIndexOfPermutation = std::numeric_limits<unsigned>::max();
  std::vector<Shapes::Vertex> bestRotation;

  const unsigned linkLimit = Shapes::rotations(shape).size();
  std::vector<unsigned> chain = {0};
  chain.reserve(32);
  std::vector<
    std::vector<Shapes::Vertex>
  > chainRotations = {identitySequence};

  /* Zero-initialize a bitset that can store if an index of permutation has
   * been discovered or not.
   */
  boost::dynamic_bitset<> discoveredIndicesOfPermutation {
    Temple::Math::factorial(Shapes::size(shape))
  };

  /* The identity sequence has been discovered (initial element of
   * chainRotations). We set bit 0 since the identity sequence has this index
   * of permutation.
   */
  discoveredIndicesOfPermutation.set(0);

  while(chain.front() < linkLimit) {
    // Generate a new rotation
    auto generatedRotation = Shapes::Properties::applyRotation(
      chainRotations.back(),
      shape,
      chain.back()
    );

    const unsigned indexOfPermutation = Temple::permutationIndex(generatedRotation);
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
      const unsigned characterIOP = Temple::permutationIndex(
        applyCharacterRotation(generatedRotation)
      );
      if(
        generatedRotation.at(reducedFusedVertex) == fusedVertex
        && characterIOP < lowestIndexOfPermutation
      ) {
        bestRotation = generatedRotation;
        lowestIndexOfPermutation = characterIOP;
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

Shapes::Vertex Composite::OrientationState::lowestEqualVertexInShape() const {
  const auto positionGroupCharacters = Shapes::Properties::positionGroupCharacters(shape);

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

  return Shapes::Vertex(findIter - std::begin(positionGroupCharacters));
}

Composite::AngleGroup Composite::OrientationState::smallestAngleGroup() const {
  // Initialize the search state
  AngleGroup angleGroup;
  angleGroup.vertices.reserve(Shapes::size(shape));
  angleGroup.angle = M_PI;

  // Go through all shape vertices excluding the fused shape position
  for(Shapes::Vertex i {0}; i < Shapes::size(shape); ++i) {
    if(i == fusedVertex) {
      continue;
    }

    const double angleToFusedPosition = Shapes::angleFunction(shape)(fusedVertex, i);

    // This naturally excludes M_PI angles from the smallest angle group
    if(fpComparator.isLessThan(angleToFusedPosition, angleGroup.angle)) {
      angleGroup.vertices = {i};
      angleGroup.angle = angleToFusedPosition;
    } else if(fpComparator.isEqual(angleToFusedPosition, angleGroup.angle) && angleGroup.angle != M_PI) {
      angleGroup.vertices.push_back(i);
    }
  }

  /* Isotropicity of rotation of the angle group can be determined as follows:
   * Group all vertices by their abstract ranking character. If, for any group
   * of vertices, the sum vector of the idealized coordinates of the shape
   * vertices does not lie along the bond line, the rotation is non-isotropic.
   */

  // Collect the vertices by ranking
  const auto rankingGroups = Temple::groupByMapping(
    angleGroup.vertices,
    Temple::Functor::at(characters)
  );

  const Eigen::ParametrizedLine<double, 3> bondLine {
    Eigen::Vector3d::Zero(),
    Shapes::coordinates(shape).col(fusedVertex)
  };

  angleGroup.isotropic = Temple::all_of(
    rankingGroups,
    [&](const std::vector<Shapes::Vertex>& vertexGroup) -> bool {
      Eigen::Vector3d vertexSum = Eigen::Vector3d::Zero();
      for(Shapes::Vertex v : vertexGroup) {
        vertexSum += Shapes::coordinates(shape).col(v);
      }

      return bondLine.squaredDistance(vertexSum) < 0.01;
    }
  );

  assert(
    std::is_sorted(
      std::begin(angleGroup.vertices),
      std::end(angleGroup.vertices)
    )
  );

  return angleGroup;
}

double Composite::perpendicularSubstituentAngle(
  const double angleFromBoundShapeVertex,
  const double angleBetweenSubstituents
) {
  assert(angleFromBoundShapeVertex != M_PI);

  return std::acos(
    1.0 - (
      1.0 - std::cos(angleBetweenSubstituents)
    ) / (
      std::pow(
        std::sin(angleFromBoundShapeVertex),
        2
      )
    )
  );
}

std::vector<Shapes::Vertex> Composite::generateRotation(
  const Shapes::Shape shape,
  const Shapes::Vertex fixedVertex,
  const std::vector<Shapes::Vertex>& changedVertices
) {
  const auto periodicities = Temple::map(
    Temple::iota<unsigned>(Shapes::rotations(shape).size()),
    [&shape](const unsigned rotationFunctionIndex) -> unsigned {
      return Shapes::Properties::rotationPeriodicity(
        shape,
        Shapes::rotations(shape).at(rotationFunctionIndex)
      );
    }
  );

  auto rotationAltersPositions = [&](const std::vector<Shapes::Vertex>& rotation) -> bool {
    return Temple::all_of(
      changedVertices,
      [&rotation](const unsigned shapeVertex) -> bool {
        return rotation.at(shapeVertex) != shapeVertex;
      }
    );
  };

  // Which rotation indicates a rotation around the bound shape vertex?
  std::vector<unsigned> rotationUses (periodicities.size(), 0);
  ++rotationUses.back();

  std::vector<Shapes::Vertex> rotation;
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

    Temple::sort(rotationIndexApplicationSequence);

    do {
      // Create the rotation using the index application sequence front-to-back
      rotation = Temple::iota<Shapes::Vertex>(Shapes::size(shape));

      for(const auto r : rotationIndexApplicationSequence) {
        rotation = Shapes::Properties::applyRotation(
          rotation,
          shape,
          r
        );
      }

      // Determine if the current rotation matches all required criteria
      if(
        rotation.at(fixedVertex) == fixedVertex
        && rotationAltersPositions(rotation)
      ) {
        rotationFound = true;
      }
    } while(
      !rotationFound
      && Temple::next_permutation(rotationIndexApplicationSequence)
    );

  } while(
    !rotationFound
    && Temple::nextCombinationPermutation(rotationUses, periodicities)
  );

  if(rotationFound) {
    return rotation;
  }

  return {};
}

std::vector<Shapes::Vertex> Composite::rotation(
  const Shapes::Shape shape,
  const Shapes::Vertex fixedVertex,
  const std::vector<Shapes::Vertex>& perpendicularPlanePositions
) {
  // Three possibilities:

  if(perpendicularPlanePositions.size() > 1) {
    /* There are multiple elements in perpendicularPlanePositions. We have to
     * generate a rotation that keeps fixedVertex fixed but rotates the
     * perpendicularPlanePositions, ideally with a periodicity equivalent to the
     * amount of shape vertices involved.
     */
    auto candidateRotation = generateRotation(
      shape,
      fixedVertex,
      perpendicularPlanePositions
    );

    // There may be multiple elements, but no rotation. Return identity
    if(candidateRotation.empty()) {
      return {Shapes::Vertex(1)};
    }

    /* Require that the periodicity of the discovered rotation is equal to the
     * number of elements being rotated. This should be a natural property of
     * generateRotation, but it's best to be sure.
     */
    assert(
      Shapes::Properties::rotationPeriodicity(shape, candidateRotation)
      == perpendicularPlanePositions.size()
    );

    return candidateRotation;
  }

  if(perpendicularPlanePositions.size() == 1) {
    /* There is a single element in perpendicularPlanePositions. The resulting
     * rotation within that shape is the identity rotation, because this
     * single index can be rotated any which way to satisfy the other side.
     */
    return {Shapes::Vertex(1)};
  }

  /* Remaining case: There are no elements in perpendicularPlanePositions. Then
   * there is no rotation, not even identity, to help in combinatorial handling
   */
  return {};
}

Composite::PerpendicularAngleGroups Composite::inGroupAngles(
  const AngleGroup& angleGroup,
  const Shapes::Shape shape
) {
  PerpendicularAngleGroups groups;

  using RecordVector = std::vector<
    std::pair<unsigned, unsigned>
  >;

  Temple::forEach(
    Temple::Adaptors::allPairs(angleGroup.vertices),
    [&](const Shapes::Vertex a, const Shapes::Vertex b) -> void {
      const double perpendicularAngle = perpendicularSubstituentAngle(
        angleGroup.angle,
        Shapes::angleFunction(shape)(a, b)
      );

      const auto findIter = std::find_if(
        std::begin(groups),
        std::end(groups),
        [&](const auto& record) -> bool {
          return Temple::any_of(
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
          RecordVector {Detail::makeOrderedPair(a, b)}
        );
      } else {
        findIter->first.push_back(perpendicularAngle);
        findIter->second.emplace_back(
          Detail::makeOrderedPair(a, b)
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
) : orientations_ { std::move(first), std::move(second) },
    alignment_(alignment)
{
  // Do not construct the ordered pair of OrientationStates with same identifier
  assert(orientations_.first.identifier != orientations_.second.identifier);

  /* In order to get meaningful indices of permutation, combinations of
   * shapes across fused positions within the same group of shape
   * vertices (e.g. equatorial or apical in square pyramidal) must be the same.
   *
   * In order to achieve this, the OrientationStates is transformed by a
   * rotation that temporarily places the fused position at the lowest index
   * shape vertex in its shape. After permutation are generated,
   * the orientation state is transformed back.
   */
  const auto firstReversionMapping = orientations_.first.transformToCanonical();
  const auto secondReversionMapping = orientations_.second.transformToCanonical();

  /* Find the group of shape vertices with the smallest angle to the
   * fused position (these are the only important ones when considering
   * relative arrangements across the bond).
   */
  const auto angleGroups = orientations_.map(
    [](const OrientationState& orientation) -> AngleGroup {
      auto angleGroup = orientation.smallestAngleGroup();
      /* Order both AngleGroups' vertices by descending ranking and
       * index to get canonical initial combinations
       */
      Temple::sort(
        angleGroup.vertices,
        [&](const unsigned a, const unsigned b) -> bool {
          return (
            std::tie(orientation.characters.at(a), a)
            > std::tie(orientation.characters.at(b), b)
          );
        }
      );
      return angleGroup;
    }
  );

  /* From the angle groups' characters, we can figure out if all
   * stereopermutations that are generated will be ranking-wise equivalent
   * spatially despite differing in the shape vertex at which the equally
   * ranked substituents are placed. This is important information for deciding
   * whether a Composite yields a stereogenic object.
   */
  isotropic_ = angleGroups.first.isotropic || angleGroups.second.isotropic;

  /* Generate a set of stereopermutations for this particular combination,
   * which can then be indexed
   *
   * Range of combinatorial possibilities:
   * - Either side has zero shape vertices in the smallest angle group:
   *   No relative positioning possible, this Composite has zero
   *   stereopermutations
   * - Both sides have one shape vertex in the smallest angle group:
   *   Dihedrals for both can be cis / trans
   * - One side has one shape vertices in the smallest angle group:
   *   Dihedral is 0 to one shape vertex of the larger side, X to the others
   * - Both sides have multiple shape vertices in the smallest angle group:
   *   Figure out the relative angles between positions in each angle group and
   *   try to find matches across groups -> these can be arranged in a coplanar
   *   fashion. Then each rotation on one side generates a new overlay
   *   possibility.
   *
   * NOTE: The central atom of both shapes is always placed at the origin
   * in the coordinate definitions.
   */
  auto firstCoordinates = Shapes::coordinates(orientations_.first.shape);
  // Rotate left fused position onto <1, 0, 0>
  Detail::rotateCoordinates(
    firstCoordinates,
    firstCoordinates.col(orientations_.first.fusedVertex).normalized(),
    Eigen::Vector3d::UnitX()
  );

  auto secondCoordinates = Shapes::coordinates(orientations_.second.shape);
  // Rotate right fused position onto <-1, 0, 0>
  Detail::rotateCoordinates(
    secondCoordinates,
    secondCoordinates.col(orientations_.second.fusedVertex).normalized(),
    -Eigen::Vector3d::UnitX()
  );

  // Translate positions by <1, 0, 0>
  Detail::translateCoordinates(
    secondCoordinates,
    Eigen::Vector3d::UnitX()
  );

  auto getDihedral = [&](const Shapes::Vertex f, const Shapes::Vertex s) -> double {
    return Detail::dihedral(
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
   * thinking up an elegant solution that can satisfy all possible shapes..
   */

  // Generate all arrangements regardless of whether the Composite is isotropic
  Temple::forEach(
    Temple::Adaptors::allPairs(
      angleGroups.first.vertices,
      angleGroups.second.vertices
    ),
    [&](const Shapes::Vertex f, const Shapes::Vertex s) -> void {
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
         * to the next shape vertex in some direction. We'll go for most
         * negative.
         */

        const auto dihedrals = Temple::map(
          angleGroups.second.vertices,
          [&](const Shapes::Vertex secondShapeVertex) -> double {
            double dihedral = getDihedral(f, secondShapeVertex);
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

      auto dihedralList = Temple::sorted(
        Temple::map(
          Temple::Adaptors::allPairs(
            angleGroups.first.vertices,
            angleGroups.second.vertices
          ),
          [&](const Shapes::Vertex a, const Shapes::Vertex b) -> DihedralTuple {
            return {a, b, getDihedral(a, b)};
          }
        )
      );

      if(
        !Temple::any_of(
          stereopermutations_,
          [&dihedralList](const auto& rhsDihedralList) -> bool {
            return Detail::dihedralClose(
              std::get<2>(dihedralList.front()),
              std::get<2>(rhsDihedralList.front()),
              Temple::Math::toRadians(15.0)
            );
          }
        )
      ) {
        stereopermutations_.push_back(std::move(dihedralList));
      }
    }
  );

  /* For situations in which only one position exists in both shapes, add
   * the trans dihedral possibility explicitly
   */
  if(
    angleGroups.first.vertices.size() == 1
    && angleGroups.second.vertices.size() == 1
  ) {
    // Add trans dihedral possibility
    stereopermutations_.emplace_back(
      std::vector<DihedralTuple> {
        DihedralTuple {
          angleGroups.first.vertices.front(),
          angleGroups.second.vertices.front(),
          M_PI
        }
      }
    );
  }

  // Revert the OrientationStates and transform the stereopermutations too
  orientations_.first.revert(firstReversionMapping);
  orientations_.second.revert(secondReversionMapping);

  for(auto& stereopermutation : stereopermutations_) {
    for(auto& dihedralTuple : stereopermutation) {
      auto findIter = std::find(
        std::begin(firstReversionMapping),
        std::end(firstReversionMapping),
        std::get<0>(dihedralTuple)
      );

      assert(findIter != std::end(firstReversionMapping));

      std::get<0>(dihedralTuple) = findIter - std::begin(firstReversionMapping);

      assert(std::get<0>(dihedralTuple) != orientations_.first.fusedVertex);

      findIter = std::find(
        std::begin(secondReversionMapping),
        std::end(secondReversionMapping),
        std::get<1>(dihedralTuple)
      );

      assert(findIter != std::end(secondReversionMapping));

      std::get<1>(dihedralTuple) = findIter - std::begin(secondReversionMapping);

      assert(std::get<1>(dihedralTuple) != orientations_.second.fusedVertex);
    }
  }

  /* Reverse the stereopermutation sequence. This is so that the indices of the
   * generated permutations yield the following simple comparison:
   *
   *   0 is E, 1 is Z
   *   1 > 0 == Z > E
   */
  std::reverse(
    std::begin(stereopermutations_),
    std::end(stereopermutations_)
  );
}

void Composite::applyIdentifierPermutation(const std::vector<std::size_t>& permutation) {
  for(auto& orientationState : orientations_) {
    orientationState.identifier = permutation.at(orientationState.identifier);
  }
}

unsigned Composite::permutations() const {
  return stereopermutations_.size();
}

Composite::Alignment Composite::alignment() const {
  return alignment_;
}

const std::vector<Composite::DihedralTuple>& Composite::dihedrals(unsigned permutationIndex) const {
  return stereopermutations_.at(permutationIndex);
}

bool Composite::isIsotropic() const {
  return isotropic_;
}

unsigned Composite::order() const {
  auto countDistinct = [&](auto&& f) {
    std::set<unsigned> positions;
    for(const DihedralTuple& t : stereopermutations_.front()) {
      positions.insert(f(t));
    }
    return positions.size();
  };

  return std::max(
    countDistinct(Temple::Functor::first),
    countDistinct(Temple::Functor::second)
  );
}

const Temple::OrderedPair<Composite::OrientationState>& Composite::orientations() const {
  return orientations_;
}

Composite::PermutationsList::const_iterator Composite::begin() const {
  return std::begin(stereopermutations_);
}

Composite::PermutationsList::const_iterator Composite::end() const {
  return std::end(stereopermutations_);
}

bool Composite::operator < (const Composite& other) const {
  return orientations_ < other.orientations_;
}

bool Composite::operator == (const Composite& other) const {
  return orientations_ == other.orientations_;
}

bool Composite::operator != (const Composite& other) const {
  return !(*this == other);
}

} // namespace Stereopermutations
} // namespace Molassembler
} // namespace Scine
