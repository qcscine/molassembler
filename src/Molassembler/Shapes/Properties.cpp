/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Shapes/Properties.h"

#include "boost/functional/hash.hpp"
#include <Eigen/Dense>
#include "Molassembler/Temple/Adaptors/Transform.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/GroupBy.h"
#include "Molassembler/Temple/Optionals.h"
#include "Molassembler/Temple/Permutations.h"
#include "Molassembler/Temple/Variadic.h"
#include "Molassembler/Temple/constexpr/Numeric.h"

#include "Molassembler/Shapes/Data.h"

#include <unordered_set>

namespace Scine {
namespace Molassembler {
namespace Shapes {
namespace Detail {

/*!
 * Generates a vector containing strictly monotonically increasing natural
 * numbers representing the range [start, end). If start == end, an empty
 * range is returned.
 */
template<typename NumericType>
std::vector<NumericType> range(
  const NumericType start,
  const NumericType end
) {
  if(start == end) {
    return {};
  }

  std::vector<NumericType> values (end - start);
  std::iota(
    values.begin(),
    values.end(),
    start
  );
  return values;
}

} // namespace Detail

namespace Properties {

std::vector<Vertex> applyPermutation(
  const Permutation& occupation,
  const Permutation& permutation
) {
  assert(occupation.size() == permutation.size());

  return Temple::map(
    permutation,
    Temple::Functor::at(occupation)
  );
}

std::vector<Vertex> applyRotation(
  const std::vector<Vertex>& occupation,
  const Shape shape,
  unsigned rotationFunctionIndex
) {
  return applyPermutation(
    occupation,
    Shapes::rotations(shape).at(rotationFunctionIndex)
  );
}

unsigned rotationPeriodicity(
  const Shapes::Shape shape,
  const std::vector<Vertex>& rotation
) {
  assert(rotation.size() == Shapes::size(shape));

  const auto initialIndices = Temple::iota<Vertex>(Vertex(Shapes::size(shape)));

  std::vector<Vertex> modified = applyPermutation(initialIndices, rotation);

  unsigned i = 1;
  for(/* */; modified != initialIndices && i < 20; ++i) {
    modified = applyPermutation(modified, rotation);
  }

  // No rotation should have a periodicity of 20.
  assert(i != 20);

  return i;
}

bool isRotation(const Permutation& rotation) {
  return Temple::sorted(rotation) == Temple::iota<Vertex>(rotation.size());
}

std::vector<std::vector<Vertex>> positionGroups(const Shapes::Shape shape) {
  /* The idea behind this algorithm is that an individual rotation (which is
   * really just an index permutation) could be interpreted as a directed graph
   * whose connected components represent sets of indices that can be
   * interchanged through multiple applications of the rotation.
   *
   * Then all that is necessary is to somehow combine multiple rotations into
   * another rotation in which overlapping connected components are merged.
   */
  const unsigned S = Shapes::size(shape);

  using IndexGroups = std::vector<Permutation>;

  // Helper function to find an index in an IndexGroup
  auto nestedFind = [](const IndexGroups& a, Vertex i) {
    return Temple::find_if(
      a,
      [&](const auto& group) -> bool {
        return Temple::find(group, i) != std::end(group);
      }
    );
  };

  // Fetches all indices that are part of a permutation's index loop starting at i
  auto loopVertices = [](
    const Permutation& rotation,
    const Vertex i
  ) -> std::vector<Vertex> {
    assert(isRotation(rotation));
    std::vector<Vertex> vertices(1, i);
    unsigned j = rotation.at(i);
    while(j != i) {
      vertices.emplace_back(j);
      j = rotation.at(j);
    }
    return vertices;
  };

  // Combine two rotations into one where overlapping connected components are merged
  auto merge = [&loopVertices](
    Permutation a,
    const Permutation& b
  ) -> Permutation {
    assert(isRotation(a) && isRotation(b));
    const unsigned R = a.size();
    assert(b.size() == R);
    for(Vertex i {0}; i < R; ++i) {
      // Self-reference is uninteresting
      if(b.at(i) == i) {
        continue;
      }

      // So now i is definitely part of a loop of more than one vertex.
      // If the permutation matches, we can already continue
      if(b.at(i) == a.at(i)) {
        continue;
      }

      // Now we have to establish which vertices are in the loop in a
      auto aLoopVertices = loopVertices(a, i);

      // If it's already a full loop, we cannot possibly merge anything else
      if(aLoopVertices.size() == R) {
        return a;
      }

      // If b.at(i) is already in that loop, we can continue.
      if(Temple::find(aLoopVertices, b.at(i)) != std::end(aLoopVertices)) {
        continue;
      }

      // Figure out which vertices from b's loop aren't yet in a
      auto bLoopVertices = Temple::sorted(loopVertices(b, i));
      Temple::sort(aLoopVertices);
      std::vector<Vertex> bVerticesNotInA;
      std::set_difference(
        std::begin(bLoopVertices),
        std::end(bLoopVertices),
        std::begin(aLoopVertices),
        std::end(aLoopVertices),
        std::back_inserter(bVerticesNotInA)
      );

      // Add in any vertices from loops in a that bVerticesNotInA have
      const unsigned t = bVerticesNotInA.size();
      for(Vertex j {0}; j < t; ++j) {
        for(Vertex k : loopVertices(a, bVerticesNotInA.at(j))) {
          if(Temple::find(bVerticesNotInA, k) == std::end(bVerticesNotInA)) {
            bVerticesNotInA.push_back(k);
          }
        }
      }

      // Found full set
      if(bVerticesNotInA.size() == R - 1) {
        a = {Vertex(R - 1)};
        for(Vertex k {0}; k < R - 1; ++k) {
          a.push_back(k);
        }
        return a;
      }

      // And insert them between i -> a.at(i)
      const unsigned hook = a.at(i);
      unsigned previous = i;
      for(unsigned j : bVerticesNotInA) {
        a.at(previous) = j;
        previous = j;
      }
      a.at(previous) = hook;

      assert(isRotation(a));
    }

    return a;
  };

  // Extract all loop groups from a rotation
  auto connectedComponents = [&](const Permutation& rotation) -> IndexGroups {
    const unsigned R = rotation.size();

    IndexGroups groups;
    for(Vertex i {0}; i < R; ++i) {
      if(nestedFind(groups, i) != std::end(groups)) {
        // Already found i
        continue;
      }

      groups.push_back(loopVertices(rotation, i));
    }

    return groups;
  };

  // Merge all rotations and then calculate the connected components
  const auto identity = Temple::iota<Vertex>(S);
  return connectedComponents(
    Temple::accumulate(rotations(shape), identity, merge)
  );
}

std::vector<char> positionGroupCharacters(const Shape shape) {
  const unsigned S = size(shape);
  // Transpose to a character-based symbolic representation
  std::vector<char> characterRepresentation (S);
  char currentChar = 'A';
  for(const auto& equalSet : positionGroups(shape)) {
    for(const auto& equalIndex : equalSet) {
      characterRepresentation.at(equalIndex) = currentChar;
    }

    ++currentChar;
  }

  return characterRepresentation;
}

std::vector<Vertex> inverseRotation(const std::vector<Vertex>& rotation) {
  const unsigned N = rotation.size();

  std::vector<Vertex> permutation (N);

  for(Vertex i {0}; i < N; ++i) {
    permutation.at(rotation.at(i)) = i;
  }

  return permutation;
}

Eigen::Vector3d getCoordinates(
  const Shapes::Shape shape,
  const boost::optional<Vertex>& vertexOption
) {
  if(vertexOption) {
    assert(vertexOption.value() < Shapes::size(shape));

    return coordinates(shape).col(
      vertexOption.value()
    );
  }

  return {0, 0, 0};
}

double getTetrahedronVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
) {
  return (i - l).dot(
    (j - l).cross(k - l)
  );
}

double calculateAngleDistortion(
  const Shapes::Shape from,
  const Shapes::Shape to,
  const std::vector<Vertex>& indexMapping
) {
  const unsigned mappingIndexLimit = std::min(
    Shapes::size(from),
    Shapes::size(to)
  );

  assert(indexMapping.size() >= mappingIndexLimit);
  assert(
    std::abs(
      static_cast<int>(Shapes::size(from)) - static_cast<int>(Shapes::size(to))
    ) <= 1
  );

  double angularDistortion = 0;

  for(Vertex i(0); i < mappingIndexLimit; ++i) {
    for(Vertex j(i + 1); j < mappingIndexLimit; ++j) {
      angularDistortion += std::fabs(
        Shapes::angleFunction(from)(i, j)
        - Shapes::angleFunction(to)(
          indexMapping.at(i),
          indexMapping.at(j)
        )
      );
    }
  }

  return angularDistortion;
}

boost::optional<Vertex> propagateIndexOptionalThroughMapping(
  const boost::optional<Vertex>& indexOptional,
  const std::vector<Vertex>& indexMapping
) {
  return Temple::Optionals::map(
    indexOptional,
    [&](const Vertex v) -> Vertex {
      return indexMapping.at(v);
    }
  );
}


double calculateChiralDistortion(
  const Shapes::Shape from,
  const Shapes::Shape to,
  const std::vector<Vertex>& indexMapping
) {

  assert(indexMapping.size() >= std::min(Shapes::size(from), Shapes::size(to)));

  double chiralDistortion = 0;

  for(const auto& tetrahedron : Shapes::tetrahedra(from)) {
    chiralDistortion += std::fabs(
      getTetrahedronVolume(
        getCoordinates(from, tetrahedron.at(0)),
        getCoordinates(from, tetrahedron.at(1)),
        getCoordinates(from, tetrahedron.at(2)),
        getCoordinates(from, tetrahedron.at(3))
      ) - getTetrahedronVolume(
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(0), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(1), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(2), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(3), indexMapping)
        )
      )
    );
  }

  return chiralDistortion;
}


std::set<
  std::vector<Vertex>
> generateAllRotations(
  const Shapes::Shape shape,
  const std::vector<Vertex>& indices
) {
  /* Replacing set here with unordered_set does not yield any speed
   * improvements. The hash function is comparatively too expensive for the
   * small number of elements in the set to set it apart from the tree.
   */
  assert(Shapes::size(shape) == indices.size());

  // Idea: Tree-like expansion of all possible combinations of rotations.
  using IndicesList = std::vector<Vertex>;

  std::set<IndicesList> allRotations = {indices};

  /* We keep a chain of all applied rotations that led to a specific index
   * sequence.  The upper limit for every link in the chain is the number of
   * rotations in the shape
   */
  unsigned linkLimit = Shapes::rotations(shape).size();

  std::vector<unsigned> chain = {0};
  /* It's also necessary to keep the structures themselves since we may
   * backtrack and apply a new rotation to a previous structure. Keeping track
   * is cheaper than applying all rotations in the chain to the initial indices
   */
  std::vector<IndicesList> chainStructures = {indices};

  /* Loop is broken when the very first link in the chain has been incremented
   * to the link limit
   */
  while(chain.front() < linkLimit) {
    // perform rotation
    // copy the last element in chainStructures
    auto generated = applyRotation(
      chainStructures.back(),
      shape,
      chain.back()
    );

    // is it something new?
    if(allRotations.count(generated) == 0) {
      // then add it to the set
      allRotations.insert(generated);

      // add it to the chain
      chainStructures.push_back(generated);
      chain.emplace_back(0);
    } else {
      // collapse the chain until we are at an incrementable position (if need be)
      while(
        chain.size() > 1 // retain at least first link in chain
        && chain.back() == linkLimit - 1 // remove link only if just below limit
      ) {
        chain.pop_back();
        chainStructures.pop_back();
      }

      // increment last position in chain
      ++chain.back();
    }
  }

  return allRotations;
}

std::vector<Vertex> applyIndexMapping(
  const Shapes::Shape to,
  const std::vector<Vertex>& mapping
) {
  /* Creates the list of indices in the target shape. Why is this necessary?
   *
   * E.g. An index mapping from linear to T-shaped. The individual
   * shape-internal vertex numbering schemes are indicated.
   *
   *  1  –▶  0
   *  |      |
   * (_)    (_) – 1 (new)
   *  |      |                Line pos. 0 to
   *  0  –▶  2                pos. 2 in Tshaped
   *                                 |  ┌– Line pos. 1 to pos. 0 in Tshaped
   *                                 |  |  ┌– This position is new
   * The mapping is represented as {2, 0, 1}.
   *
   * This function writes the indices of the original mapping into the target
   * shape's indexing scheme.
   *
   * For this example, this returns {1, 2, 0}:
   *
   *  1 (at pos 0 in internal indexing scheme)
   *  |
   * (_) – 2 (etc.)
   *  |
   *  0
   *
   * For this example, this returns {1, 2, 0}.
   *
   * The closely related mapping {0, 2, 1} yields target indices {0, 2, 1}.
   *
   * Which of these properties are related by target shape rotations?
   *
   *
   *     mapping       target indices
   * -----------------------------------
   *    {2, 0, 1}   =>   {1, 2, 0}
   *        ▲                ▲
   *        |                |
   *        X                | C2 rotation in T-shape shape
   *        |                |
   *        ▼                ▼
   *    {0, 2, 1}   =>   {0, 2, 1}
   *
   */
  std::vector<Vertex> shapeVertices (Shapes::size(to));

  for(Vertex i {0}; i < Shapes::size(to); ++i) {
    shapeVertices.at(
      mapping.at(i)
    ) = i;
  }

  return shapeVertices;
}

DistortionInfo::DistortionInfo(
  std::vector<Vertex> passIndexMapping,
  const double passAngularDistortion,
  const double passChiralDistortion
) : indexMapping(std::move(passIndexMapping)),
    angularDistortion(passAngularDistortion),
    chiralDistortion(passChiralDistortion)
{}

std::size_t hash_value(const std::vector<unsigned>& permutation) {
  return Temple::permutationIndex(permutation);
}

std::vector<DistortionInfo> shapeTransitionMappings(
  const Shape from,
  const Shape to
) {

  /* Symmetries must be adjacent in size (0 = rearrangement,
   * +1 = ligand gain. Ligand loss is a special case where a specific position
   * in the shape group is removed, and is not covered here!
   */
  assert(
    (std::set<int> {0, 1}).count(
      static_cast<int>(Shapes::size(to))
      - static_cast<int>(Shapes::size(from))
    ) == 1
  );

  /* Base idea: We need to go through all possible mappings. In situations where
   * the target shape has one more or one fewer ligand, the last index in
   * the current sequence is either the added or removed ligand, and merely the
   * others are used to calculate the angular and chiral distortions involved
   * in the transition.
   *
   * For instance, from linear to T-shaped
   *
   *   0 - ( ) - 1     ->   0 - ( ) - 2
   *                             |
   *                             1
   *
   *   A mapping of {0 1 2} means that the new ligand is inserted at the 2
   *   position of T-shaped, which would involve distorting the 0-1 angle by
   *   90°. The optimal mapping would be {0 2 1} (or its equivalent rotation
   *   {2 0 1}), which does not have any angular distortion.
   */

  const unsigned largerSize = std::max(
    Shapes::size(from),
    Shapes::size(to)
  );

  auto indexMapping = Temple::iota<Vertex>(largerSize);

  // Store all distortions calculated
  std::vector<DistortionInfo> distortions;

  // Need to keep track of rotations of mappings to avoid repetition
  using PermutationType = std::vector<Vertex>;
  std::unordered_set<
    PermutationType,
    boost::hash<PermutationType>
  > encounteredShapeMappings;

  /* Using std::next_permutation generates all possible mappings!
   * do-while is required since the very first mapping must be considered before
   * calling next_permutation for the first time
   */
  do {
    if( // is the mapping new?
      encounteredShapeMappings.count(
        applyIndexMapping(to, indexMapping)
      ) == 0
    ) {
      /* Add it to the set of possible distortions, calculate angular and
       * chiral distortion involved in the mapping
       */
      distortions.emplace_back(
        indexMapping,
        calculateAngleDistortion(from, to, indexMapping),
        calculateChiralDistortion(from, to, indexMapping)
      );

      /* Any rotations of the mapping in the target shape are equivalent, we
       * do not want to count these as an additional multiplicity, so we
       * generate them and add them to the encountered mappings
       */
      auto allRotations = generateAllRotations(
        to,
        applyIndexMapping(to, indexMapping)
      );

      encounteredShapeMappings.insert(
        std::make_move_iterator(allRotations.begin()),
        std::make_move_iterator(allRotations.end())
      );
    }
  } while (std::next_permutation(indexMapping.begin(), indexMapping.end()));

  return distortions;
}

std::vector<DistortionInfo> ligandLossTransitionMappings(
  const Shapes::Shape from,
  const Shapes::Shape to,
  const Vertex positionInSourceShape
) {
  // Ensure we are dealing with ligand loss
  assert(Shapes::size(to) + 1 == Shapes::size(from));
  assert(positionInSourceShape < Shapes::size(from));

  /* Generate the index mapping specific to this position loss in the target
   * shape.
   *
   * The possible distortions are determined from the equivalent case of adding
   * the ligand (which was to be deleted) to the smaller shape.
   *
   * The deleted index is added to the end. The lowest permutation of the mapping
   * from the smaller to the larger shape is then:
   *
   * 0, 1, ..., (pos - 1), (pos + 1), ..., (from_size - 1), pos
   *
   * NOTE: pos simply isn't added to the mapping (the distortion component
   * calculating functions tolerate this), simplifying the permutation
   */
  std::vector<Vertex> indexMapping = Temple::variadic::concatenate(
    Temple::iota<Vertex>(positionInSourceShape),
    Detail::range(Vertex(positionInSourceShape + 1), Vertex(Shapes::size(from)))
  );

  /* NOTE: From here the algorithm is identical to shapeTransitionMappings
   * save that to and from are swapped in all occasions
   * and that std::next_permutation is only called on the subset excluding the
   * last position (the one that is added / deleted).
   *
   * AAAND there's no applyTransitionMapping for some reason
   */

  std::vector<DistortionInfo> distortions;

  std::set<
    std::vector<Vertex>
  > encounteredShapeMappings;

  do {
    if(encounteredShapeMappings.count(indexMapping) == 0) {
      distortions.emplace_back(
        indexMapping,
        calculateAngleDistortion(to, from, indexMapping),
        calculateChiralDistortion(to, from, indexMapping)
      );

      auto allRotations = generateAllRotations(to, indexMapping);

      encounteredShapeMappings.insert(
        std::make_move_iterator(allRotations.begin()),
        std::make_move_iterator(allRotations.end())
      );
    }
  } while (std::next_permutation(indexMapping.begin(), indexMapping.end()));

  return distortions;
}

ShapeTransitionGroup::ShapeTransitionGroup(
  std::vector<
    std::vector<Vertex>
  > passIndexMappings,
  const double passAngleDistortion,
  const double passChiralDistortion
) : indexMappings(std::move(passIndexMappings)),
    angularDistortion(passAngleDistortion),
    chiralDistortion(passChiralDistortion)
{}

ShapeTransitionGroup selectBestTransitionMappings(
  const std::vector<DistortionInfo>& distortions
) {
  /* We are interested only in the those transitions that have the very lowest
   * angular distortion, and within that set only the lowest chiral distortion,
   * so we sub-select within the generated set
   */

  const double lowestAngularDistortion = std::min_element(
    std::begin(distortions),
    std::end(distortions),
    [](const auto& a, const auto& b) -> bool {
      return a.angularDistortion < b.angularDistortion;
    }
  )->angularDistortion;

  std::vector<unsigned> viableDistortions;
  for(unsigned i = 0; i < distortions.size(); ++i) {
    if(distortions.at(i).angularDistortion < lowestAngularDistortion + floatingPointEqualityThreshold) {
      viableDistortions.push_back(i);
    }
  }

  // And now sub-set further on the lowest chiral distortion
  const double lowestChiralDistortion = distortions.at(
    *std::min_element(
      std::begin(viableDistortions),
      std::end(viableDistortions),
      [&](const unsigned a, const unsigned b) -> bool {
        return distortions.at(a).chiralDistortion < distortions.at(b).chiralDistortion;
      }
    )
  ).chiralDistortion;

  Temple::remove_if(
    viableDistortions,
    [&](const unsigned a) -> bool {
      return distortions.at(a).chiralDistortion > lowestChiralDistortion + floatingPointEqualityThreshold;
    }
  );

  return ShapeTransitionGroup(
    Temple::map(viableDistortions, [&](auto i) { return distortions.at(i).indexMapping; }),
    lowestAngularDistortion,
    lowestChiralDistortion
  );
}

unsigned numUnlinkedStereopermutations(
  const Shapes::Shape shape,
  const unsigned nIdenticalLigands
) {
  unsigned count = 1;

  auto indices = Temple::iota<Vertex>(Shapes::size(shape));

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  auto initialRotations = generateAllRotations(shape, indices);

  std::set<decltype(indices)> rotations {
    std::make_move_iterator(std::begin(initialRotations)),
    std::make_move_iterator(std::end(initialRotations))
  };

  while(std::next_permutation(indices.begin(), indices.end())) {
    if(rotations.count(indices) == 0) {
      auto allRotations = generateAllRotations(shape, indices);
      for(const auto& rotation : allRotations) {
        rotations.insert(rotation);
      }

      ++count;
    }
  }

  return count;
}

bool hasMultipleUnlinkedStereopermutations(
  const Shapes::Shape shape,
  const unsigned nIdenticalLigands
) {
  if(nIdenticalLigands == Shapes::size(shape)) {
    return false;
  }

  auto indices = Temple::iota<Vertex>(Shapes::size(shape));

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  auto initialRotations = generateAllRotations(shape, indices);

  std::set<decltype(indices)> rotations {
    std::make_move_iterator(std::begin(initialRotations)),
    std::make_move_iterator(std::end(initialRotations))
  };

  while(Temple::next_permutation(indices)) {
    if(rotations.count(indices) == 0) {
      return true;
    }
  }

  return false;
}

Shape mostSymmetric(std::vector<Shape> selection) {
  std::sort(
    std::begin(selection),
    std::end(selection),
    [](const Shape a, const Shape b) -> bool {
      return std::make_tuple(
        rotations(a).size(),
        nameIndex(b)
      ) < std::make_tuple(
        rotations(b).size(),
        nameIndex(a)
      );
    }
  );

  return selection.back();
}

Shape mostSymmetric(const unsigned shapeSize) {
  std::vector<Shape> propositions;
  propositions.reserve(8);

  for(const Shape proposition : allShapes) {
    if(size(proposition) == shapeSize) {
      propositions.push_back(proposition);
    }
  }

  return mostSymmetric(std::move(propositions));
}

} // namespace Properties
} // namespace Shapes
} // namespace Molassembler
} // namespace Scine
