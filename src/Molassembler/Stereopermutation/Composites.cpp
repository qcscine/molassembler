/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutation/Composites.h"

#include "boost/dynamic_bitset.hpp"
#include "Eigen/Geometry"

#include "Molassembler/Shapes/Properties.h"
#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Detail/Cartesian.h"
#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Adaptors/CyclicFrame.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Functor.h"
#include "Molassembler/Temple/GroupBy.h"
#include "Molassembler/Temple/Permutations.h"
#include "Molassembler/Temple/Stl17.h"

namespace Scine {
namespace Molassembler {
namespace Stereopermutations {
namespace {

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

void rotateCoordinates(
  Eigen::Ref<Shapes::Coordinates> positions,
  const Eigen::Vector3d& axis,
  const double angle
) {
  const unsigned cols = positions.cols();
  for(unsigned i = 0; i < cols; ++i) {
    positions.col(i) = Eigen::AngleAxisd(angle, axis) * positions.col(i);
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

} // namespace

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

bool Composite::Permutation::close(const std::vector<DihedralTuple>& comparisonDihedrals) const {
  return Temple::all_of(
    comparisonDihedrals,
    [&](const DihedralTuple& dihedral) -> bool {
      return Temple::find_if(
        dihedrals,
        [&](const DihedralTuple& existingDihedral) -> bool {
          return (
            std::get<0>(dihedral) == std::get<0>(existingDihedral)
            && std::get<1>(dihedral) == std::get<1>(existingDihedral)
            && Cartesian::dihedralDifference(std::get<2>(dihedral), std::get<2>(existingDihedral)) < 1e-8
          );
        }
      ) != std::end(dihedrals);
    }
  );
}

Composite::PermutationGenerator::PermutationGenerator(
  Temple::OrderedPair<OrientationState> passOrientations
) : orientations(std::move(passOrientations)) {
  /* In order to get meaningful indices of permutation, combinations of
   * shapes across fused positions within the same group of shape
   * vertices (e.g. equatorial or apical in square pyramidal) must be the same.
   *
   * In order to achieve this, the OrientationStates is transformed by a
   * rotation that temporarily places the fused position at the lowest index
   * shape vertex in its shape. After permutation are generated,
   * the orientation state is transformed back.
   */
  reversionMappings = std::make_pair(
    orientations.first.transformToCanonical(),
    orientations.second.transformToCanonical()
  );

  /* Find the group of shape vertices with the smallest angle to the
   * fused position (these are the only important ones when considering
   * relative arrangements across the bond).
   */
  angleGroups = orientations.map(
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
  coordinates = orientations.map(
    [](const OrientationState& orientation) {
      return Shapes::coordinates(orientation.shape);
    }
  );

  // Rotate left fused position onto <1, 0, 0>
  rotateCoordinates(
    coordinates.first,
    coordinates.first.col(orientations.first.fusedVertex).normalized(),
    Eigen::Vector3d::UnitX()
  );

  // Rotate right fused position onto <-1, 0, 0>
  rotateCoordinates(
    coordinates.second,
    coordinates.second.col(orientations.second.fusedVertex).normalized(),
    -Eigen::Vector3d::UnitX()
  );

  // Translate right positions by <1, 0, 0>
  translateCoordinates(
    coordinates.second,
    Eigen::Vector3d::UnitX()
  );
}

Composite::PermutationsList
Composite::PermutationGenerator::generateEclipsedOrStaggered(
  const Alignment alignment,
  const double deduplicationDegrees
) {
  assert(
    alignment == Alignment::Eclipsed
    || alignment == Alignment::Staggered
  );

  PermutationsList stereopermutations;

  auto findRankingEquivalent = [&](const Permutation& permutation) -> boost::optional<Permutation::VertexPair> {
    auto replaceVertices = [&](Shapes::Vertex i, Shapes::Vertex j, double dihedral) {
      return std::make_tuple(
        orientations.first.characters.at(i),
        orientations.second.characters.at(j),
        dihedral
      );
    };

    const auto characterDihedrals = Temple::map(permutation.dihedrals, replaceVertices);

    const auto equivalentIter = Temple::find_if(
      stereopermutations,
      [&](const Permutation& searchPermutation) -> bool {
        if(searchPermutation.rankingEquivalentTo) {
          return false;
        }

        if(searchPermutation.alignment != permutation.alignment) {
          return false;
        }

        return Temple::all_of(
          Temple::map(searchPermutation.dihedrals, replaceVertices),
          [&](const auto& characterTuple) -> bool {
            return Temple::find_if(
              characterDihedrals,
              [&](const auto& baseTuple) -> bool {
                return (
                  std::get<0>(characterTuple) == std::get<0>(baseTuple)
                  && std::get<1>(characterTuple) == std::get<1>(baseTuple)
                  && Cartesian::dihedralDifference(std::get<2>(characterTuple), std::get<2>(baseTuple)) <= 1e-5
                );
              }
            ) != std::end(characterDihedrals);
          }
        );
      }
    );

    if(equivalentIter != std::end(stereopermutations)) {
      return equivalentIter->alignedVertices;
    }

    return boost::none;
  };

  /* Sequentially align every pair.
   *
   * This is essentially brute-forcing the problem. I'm having a hard time
   * thinking up an elegant solution that can satisfy all possible shapes...
   */
  Temple::forEach(
    Temple::Adaptors::allPairs(
      angleGroups.first.vertices,
      angleGroups.second.vertices
    ),
    [&](const Shapes::Vertex f, const Shapes::Vertex s) -> void {
      auto permutation = align(f, s, alignment);
      if(!isDuplicate(permutation, stereopermutations, deduplicationDegrees)) {
        permutation.rankingEquivalentTo = findRankingEquivalent(permutation);
        stereopermutations.push_back(std::move(permutation));
      }
    }
  );

  /* For situations in which only one position exists in both shapes, add
   * the trans dihedral possibility explicitly. But what dihedral is set as
   * the cis dihedral dependends on the alignment.
   */
  if(
    angleGroups.first.vertices.size() == 1
    && angleGroups.second.vertices.size() == 1
  ) {
    // Add trans dihedral possibility
    stereopermutations.push_back(stereopermutations.front());
    if(alignment == Alignment::Eclipsed) {
      std::get<2>(stereopermutations.back().dihedrals.front()) = M_PI;
    } else {
      std::get<2>(stereopermutations.back().dihedrals.front()) = 0.0;
    }

    /* In staggered alignments of two single vertex sides, staggering with
     * alignment yields a pi dihedral, and the code before a zero dihedral.
     * But staggered here means offset by pi/2:
     */
    if(alignment == Alignment::Staggered) {
      for(auto& permutation : stereopermutations) {
        for(auto& dihedral : permutation.dihedrals) {
          std::get<2>(dihedral) = Cartesian::signedDihedralAngle(std::get<2>(dihedral) + M_PI / 2);
        }
      }
    }
  }

  return stereopermutations;
}

Composite::Permutation Composite::PermutationGenerator::align(
  const Shapes::Vertex firstVertex,
  const Shapes::Vertex secondVertex,
  const Alignment alignment
) {
  const Eigen::Vector3d xAxis = Eigen::Vector3d::UnitX();
  const double alignAngle = dihedral(firstVertex, secondVertex);

  // Twist the right coordinates around x so that f is ecliptic with r
  rotateCoordinates(coordinates.second, xAxis, -alignAngle);

  // Make sure the rotation leads to a zero dihedral
  assert(std::fabs(dihedral(firstVertex, secondVertex)) < 1e-10);

  // Offset the dihedral angle for staggered alignments
  if(alignment == Alignment::Staggered) {
    /* The offset angle for a staggered arrangement is half of the angle
     * to the next shape vertex in some direction. We'll go for nearest
     * negative dihedral.
     */
    const double nearestNegativeDihedral = Temple::accumulate(
      angleGroups.second.vertices,
      std::numeric_limits<double>::lowest(),
      [&](const double carry, const Shapes::Vertex secondShapeVertex) -> double {
        double angle = dihedral(firstVertex, secondShapeVertex);
        if(angle >= -1e-10) {
          angle -= 2 * M_PI;
        }

        return std::max(carry, angle);
      }
    );

    rotateCoordinates(coordinates.second, xAxis, nearestNegativeDihedral / 2);
  }

  auto vertices = std::make_pair(firstVertex, secondVertex);
  auto sortedDihedrals = Temple::sorted(
    Temple::map(
      Temple::Adaptors::allPairs(
        angleGroups.first.vertices,
        angleGroups.second.vertices
      ),
      [&](const Shapes::Vertex a, const Shapes::Vertex b) -> Permutation::DihedralTuple {
        return {a, b, dihedral(a, b)};
      }
    )
  );

  return Permutation {
    std::move(vertices),
    alignment,
    std::move(sortedDihedrals),
    boost::none
  };
}

bool Composite::PermutationGenerator::isDuplicate(
  Permutation permutation,
  const PermutationsList& permutations,
  const double degrees
) {
  return Temple::any_of(
    permutations,
    [&](const auto& existingPermutation) -> bool {
      return Cartesian::dihedralDifference(
        std::get<2>(permutation.dihedrals.front()),
        std::get<2>(existingPermutation.dihedrals.front())
      ) <= Temple::Math::toRadians(degrees);
    }
  );
}

Composite::PermutationsList Composite::PermutationGenerator::deduplicate(PermutationsList permutations, const double degrees) {
  PermutationsList deduplicated;
  for(auto& permutation : permutations) {
    if(!isDuplicate(permutation, deduplicated, degrees)) {
      deduplicated.push_back(std::move(permutation));
    }
  }
  return deduplicated;
}

Composite::PermutationsList Composite::PermutationGenerator::generate(
  const Alignment alignment,
  const double deduplicationDegrees
) {
  PermutationsList stereopermutations;

  if(alignment == Alignment::Eclipsed || alignment == Alignment::Staggered) {
    stereopermutations = generateEclipsedOrStaggered(alignment, deduplicationDegrees);
  } else {
    /* Other cases are EclipsedAndStaggered and BetweenEclipsedAndStaggered,
     * where we need both!
     */
    stereopermutations = generateEclipsedOrStaggered(Alignment::Eclipsed, deduplicationDegrees);
    auto staggeredStereopermutations = generateEclipsedOrStaggered(Alignment::Staggered, deduplicationDegrees);

    std::copy_if(
      std::make_move_iterator(std::begin(staggeredStereopermutations)),
      std::make_move_iterator(std::end(staggeredStereopermutations)),
      std::back_inserter(stereopermutations),
      [&](const auto& stereopermutation) -> bool {
        return !isDuplicate(stereopermutation, stereopermutations, deduplicationDegrees);
      }
    );

    // Sort stereopermutations starting at the dominant dihedral zero
    Temple::sort(
      stereopermutations,
      [](const Permutation& lhs, const Permutation& rhs) {
        return (
          Cartesian::positiveDihedralAngle(std::get<2>(lhs.dihedrals.front()))
          < Cartesian::positiveDihedralAngle(std::get<2>(rhs.dihedrals.front()))
        );
      }
    );

    if(alignment == Alignment::BetweenEclipsedAndStaggered) {
      // Average the dihedral angle of successive pairs
      stereopermutations = Temple::map(
        Temple::Adaptors::cyclicFrame<2>(stereopermutations),
        [](const Permutation& lhs, const Permutation& rhs) {
          auto dihedrals = Temple::map(
            Temple::Adaptors::zip(lhs.dihedrals, rhs.dihedrals),
            [](const Permutation::DihedralTuple& a, const Permutation::DihedralTuple& b) {
              assert(std::get<0>(a) == std::get<0>(b));
              assert(std::get<1>(a) == std::get<1>(b));
              return Permutation::DihedralTuple {
                std::get<0>(a),
                std::get<1>(a),
                Cartesian::dihedralAverage(std::get<2>(a), std::get<2>(b))
              };
            }
          );

          return Permutation {
            lhs.alignedVertices,
            lhs.alignment,
            std::move(dihedrals),
            lhs.rankingEquivalentTo
          };
        }
      );
    }
  }

  // Transform the stereopermutations back through the reversion mappings
  const auto revert = [](const Shapes::Vertex v, const Shapes::Permutation& mapping) {
    const auto findIter = Temple::find(mapping, v);
    assert(findIter != std::end(mapping));
    unsigned index = findIter - std::begin(mapping);
    return Shapes::Vertex {index};
  };

  for(auto& stereopermutation : stereopermutations) {
    stereopermutation.alignedVertices = std::make_pair(
      revert(stereopermutation.alignedVertices.first, reversionMappings.first),
      revert(stereopermutation.alignedVertices.second, reversionMappings.second)
    );
    for(auto& dihedralTuple : stereopermutation.dihedrals) {
      std::get<0>(dihedralTuple) = revert(std::get<0>(dihedralTuple), reversionMappings.first);
      std::get<1>(dihedralTuple) = revert(std::get<1>(dihedralTuple), reversionMappings.second);
    }
    if(stereopermutation.rankingEquivalentTo) {
      stereopermutation.rankingEquivalentTo = std::make_pair(
        revert(stereopermutation.rankingEquivalentTo->first, reversionMappings.first),
        revert(stereopermutation.rankingEquivalentTo->second, reversionMappings.second)
      );
    }
  }

  /* Ensure that each permutation is either a ranking-unique, or its referenced
   * base exists
   */
  assert(
    Temple::all_of(
      stereopermutations,
      [&](const Permutation& permutation) {
        if(!permutation.rankingEquivalentTo) {
          return true;
        }

        return Temple::find_if(
          stereopermutations,
          [&](const Permutation& searchPerm) {
            return searchPerm.alignedVertices == permutation.rankingEquivalentTo.value();
          }
        ) != std::end(stereopermutations);
      }
    )
  );

  return stereopermutations;
}

double Composite::PermutationGenerator::dihedral(
  const Shapes::Vertex firstVertex,
  const Shapes::Vertex secondVertex
) const {
  return Cartesian::dihedral(
    coordinates.first.col(firstVertex),
    Eigen::Vector3d::Zero(),
    Eigen::Vector3d::UnitX(),
    coordinates.second.col(secondVertex)
  );
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
      std::pow(std::sin(angleFromBoundShapeVertex), 2)
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

Composite::Composite(
  OrientationState first,
  OrientationState second,
  const Alignment alignment
) : orientations_ { std::move(first), std::move(second) },
    alignment_(alignment)
{
  // Do not construct the ordered pair of OrientationStates with same identifier
  assert(orientations_.first.identifier != orientations_.second.identifier);

  PermutationGenerator generator(orientations_);
  stereopermutations_ = generator.generate(alignment);

  if(alignment == Alignment::Eclipsed) {
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
}

void Composite::applyIdentifierPermutation(const std::vector<std::size_t>& permutation) {
  for(auto& orientationState : orientations_) {
    orientationState.identifier = permutation.at(orientationState.identifier);
  }
}

const Composite::PermutationsList& Composite::allPermutations() const {
  return stereopermutations_;
}

Composite::Alignment Composite::alignment() const {
  return alignment_;
}

unsigned Composite::rankingEquivalentBase(const unsigned permutation) const {
  const Permutation& stereopermutation = stereopermutations_.at(permutation);

  if(!stereopermutation.rankingEquivalentTo) {
    return permutation;
  }

  auto findIter = Temple::find_if(
    stereopermutations_,
    [&](const auto& searchPermutation) {
      return (
        searchPermutation.alignment == stereopermutation.alignment
        && searchPermutation.alignedVertices == stereopermutation.rankingEquivalentTo.value()
      );
    }
  );

  assert(findIter != std::end(stereopermutations_));

  return findIter - std::begin(stereopermutations_);
}

std::vector<unsigned> Composite::nonEquivalentPermutationIndices() const {
  std::vector<unsigned> indices;
  const unsigned N = stereopermutations_.size();
  for(unsigned i = 0; i < N; ++i) {
    if(!stereopermutations_.at(i).rankingEquivalentTo) {
      indices.push_back(i);
    }
  }
  return indices;
}

unsigned Composite::countNonEquivalentPermutations() const {
  return Temple::accumulate(
    stereopermutations_,
    0u,
    [](unsigned carry, const Permutation& permutation) -> unsigned {
      if(permutation.rankingEquivalentTo) {
        return carry;
      }

      return carry + 1;
    }
  );
}

bool Composite::isIsotropic() const {
  return countNonEquivalentPermutations() <= 1;
}

std::pair<unsigned, unsigned> Composite::orders() const {
  /* Angle groups information is gone, but we can deduce the information
   * by counting the distinct vertices at both sides of the dihedral lists of a
   * permutation.
   */
  const auto countDistinct = [&](auto&& f) {
    std::set<unsigned> positions;
    for(const Permutation::DihedralTuple& t : stereopermutations_.front().dihedrals) {
      positions.insert(f(t));
    }
    return positions.size();
  };

  return std::make_pair(
    countDistinct(Temple::Functor::first),
    countDistinct(Temple::Functor::second)
  );
}

unsigned Composite::order() const {
  const auto both = orders();
  return std::max(both.first, both.second);
}


unsigned Composite::rotationalAxisSymmetryOrder() const {
  /* NOTE: Maybe this can be replaced with the integer division
   * stereopermutations_.count() / countNonEquivalentPermutations()?
   *
   * (well, with a special case for zeroes)
   */

  std::set<unsigned> counts;
  for(const Permutation& permutation : stereopermutations_) {
    if(permutation.rankingEquivalentTo) {
      continue;
    }

    auto count = Temple::accumulate(
      stereopermutations_,
      1u,
      [&](unsigned carry, const Permutation& other) -> unsigned {
        if(
          permutation.alignment == other.alignment
          && permutation.alignedVertices == other.rankingEquivalentTo
        ) {
          return carry + 1;
        }

        return carry;
      }
    );

    counts.insert(count);
  }

  // If there are different counts, then this cannot have a rotational symmetry
  if(counts.empty()) {
    return 0;
  }

  if(counts.size() > 1) {
    return 1;
  }

  return *std::begin(counts);
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
