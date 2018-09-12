#include "chemical_symmetries/DynamicProperties.h"
#include "chemical_symmetries/Symmetries.h"

#include <Eigen/Dense>

#include "temple/Adaptors/Transform.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Functional.h"
#include "temple/GroupBy.h"
#include "temple/Variadic.h"
#include "temple/VectorView.h"

/* TODO
 * - Consider unordered_set for rotation contains, using the hash function from
 *   constexpr work
 */

namespace Symmetry {

namespace detail {

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

} // namespace detail

namespace properties {

std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const std::vector<unsigned>& rotation
) {
  assert(indices.size() == rotation.size());

  std::vector<unsigned> retv;
  retv.reserve(indices.size());

  for(const auto& index : rotation) {
    retv.push_back(
      indices.at(index)
    );
  }

  return retv;
}

std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const Symmetry::Name symmetryName,
  unsigned rotationFunctionIndex
) {
  std::vector<unsigned> retv;

  for(
    const auto& index :
    Symmetry::rotations(symmetryName).at(rotationFunctionIndex)
  ) {
    retv.push_back(
      indices.at(index)
    );
  }

  return retv;
}

unsigned rotationPeriodicity(
  const Symmetry::Name symmetryName,
  const std::vector<unsigned>& rotation
) {
  assert(rotation.size() == Symmetry::size(symmetryName));

  const auto initialIndices = temple::iota<unsigned>(Symmetry::size(symmetryName));

  std::vector<unsigned> modified = applyRotation(initialIndices, rotation);

  unsigned i = 1;
  for(/* */; modified != initialIndices && i < 20; ++i) {
    modified = applyRotation(modified, rotation);
  }

  // No rotation should have a periodicity of 20.
  assert(i != 20);

  return i;
}

std::vector<char> positionGroups(const Symmetry::Name symmetryName) {
  const unsigned S = Symmetry::size(symmetryName);

  std::vector<
    std::vector<double>
  > allAngles (S);

  // Generate sorted lists of all cross-angles for each starting index
  for(unsigned i = 0; i < S; ++i) {
    allAngles.at(i).resize(S);
    for(unsigned j = 0; j < S; ++j) {
      allAngles.at(i).at(j) = Symmetry::angleFunction(symmetryName)(i, j);
    }

    std::sort(
      std::begin(allAngles.at(i)),
      std::end(allAngles.at(i))
    );
  }

  // Group the individual symmetry indices by identical cross-angle sets
  auto groups = temple::groupByEquality(
    temple::iota<unsigned>(S),
    [&allAngles](const unsigned i, const unsigned j) -> bool {
      return allAngles.at(i) == allAngles.at(j);
    }
  );

  // Transpose to a simple symbolic representation
  std::vector<char> characterRepresentation (S);
  char currentChar = 'A';
  for(const auto& equalSet : groups) {
    for(const auto& equalIndex : equalSet) {
      characterRepresentation.at(equalIndex) = currentChar;
    }

    ++currentChar;
  }

  return characterRepresentation;
}

std::vector<unsigned> inverseRotation(const std::vector<unsigned>& rotation) {
  const unsigned N = rotation.size();

  std::vector<unsigned> permutation (N);

  for(unsigned i = 0; i < N; ++i) {
    permutation.at(rotation.at(i)) = i;
  }

  return permutation;
}

std::vector<unsigned> invertedSequence(const Symmetry::Name symmetryName) {
  auto occupation = temple::iota<unsigned>(Symmetry::size(symmetryName));

  auto getCoordinates = [](
    const Symmetry::Name symmetry,
    const boost::optional<unsigned>& indexOption
  ) -> Eigen::Vector3d {
    if(indexOption) {
      return symmetryData().at(symmetry).coordinates.at(indexOption.value());
    }

    return {0, 0, 0};
  };

  auto mapping = [&occupation](const boost::optional<unsigned>& indexOption) -> boost::optional<unsigned> {
    if(indexOption) {
      return occupation.at(indexOption.value());
    }

    return boost::none;
  };

  // Brute-force approach
  while(
    !temple::all_of(
      Symmetry::tetrahedra(symmetryName),
      [&](const auto& tetrahedronDefinition) -> bool {
        return getTetrahedronVolume(
          getCoordinates(symmetryName, mapping(tetrahedronDefinition[0])),
          getCoordinates(symmetryName, mapping(tetrahedronDefinition[1])),
          getCoordinates(symmetryName, mapping(tetrahedronDefinition[2])),
          getCoordinates(symmetryName, mapping(tetrahedronDefinition[3]))
        ) < 0;
      }
    ) && std::next_permutation(
      std::begin(occupation),
      std::end(occupation)
    )
  ) {}

  return occupation;
}

Eigen::Vector3d getCoordinates(
  const Symmetry::Name symmetryName,
  const boost::optional<unsigned>& indexInSymmetryOption
) {
  assert(
    (
      indexInSymmetryOption
      && indexInSymmetryOption.value() < Symmetry::size(symmetryName)
    ) || !indexInSymmetryOption
  );

  if(indexInSymmetryOption) {
    return Symmetry::symmetryData().at(symmetryName).coordinates.at(
      indexInSymmetryOption.value()
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
  const Symmetry::Name from,
  const Symmetry::Name to,
  const std::vector<unsigned>& indexMapping
) {
  const unsigned mappingIndexLimit = std::min(
    Symmetry::size(from),
    Symmetry::size(to)
  );

  assert(indexMapping.size() >= mappingIndexLimit);
  assert(
    std::abs(
      static_cast<int>(Symmetry::size(from)) - static_cast<int>(Symmetry::size(to))
    ) <= 1
  );

  double angularDistortion = 0;

  for(unsigned i = 0; i < mappingIndexLimit; ++i) {
    for(unsigned j = i + 1; j < mappingIndexLimit; ++j) {
      angularDistortion += std::fabs(
        Symmetry::angleFunction(from)(i, j)
        - Symmetry::angleFunction(to)(
          indexMapping.at(i),
          indexMapping.at(j)
        )
      );
    }
  }

  return angularDistortion;
}

unsigned long hashIndexList(const std::vector<unsigned>& indexList) {
  constexpr unsigned maxDigitsStoreable = temple::Math::floor(
    temple::Math::log10(
      static_cast<double>(
        std::numeric_limits<unsigned long>::max()
      )
    )
  );

  if(indexList.size() > maxDigitsStoreable) {
    throw std::logic_error("hashIndexList cannot hash such a long index list");
  }

  long unsigned hash = 0;
  unsigned tenPowers = 1;

  for(const auto index: indexList) {
    if(index > 9) {
      throw std::logic_error("hashIndexList: Index list contains numbers greater than 9!");
    }
    hash += tenPowers * index;
    tenPowers *= 10;
  }

  return hash;
}

boost::optional<unsigned> propagateIndexOptionalThroughMapping(
  const boost::optional<unsigned>& indexOptional,
  const std::vector<unsigned>& indexMapping
) {
  if(indexOptional) {
    return indexMapping.at(indexOptional.value());
  }

  return boost::none;
}


double calculateChiralDistortion(
  const Symmetry::Name from,
  const Symmetry::Name to,
  const std::vector<unsigned>& indexMapping
) {

  assert(
    indexMapping.size() >= std::min(
      Symmetry::size(from),
      Symmetry::size(to)
    )
  );

  double chiralDistortion = 0;

  for(const auto& tetrahedron : Symmetry::tetrahedra(from)) {
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
  std::vector<unsigned>
> generateAllRotations(
  const Symmetry::Name symmetryName,
  const std::vector<unsigned>& indices
) {
  assert(Symmetry::size(symmetryName) == indices.size());

  // Idea: Tree-like expansion of all possible combinations of rotations.
  using IndicesList = std::vector<unsigned>;

  std::set<IndicesList> allRotations = {indices};

  /* We keep a chain of all applied rotations that led to a specific index
   * sequence.  The upper limit for every link in the chain is the number of
   * rotations in the symmetry
   */
  unsigned linkLimit = Symmetry::rotations(symmetryName).size();

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
      symmetryName,
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

std::vector<unsigned> applyIndexMapping(
  const Symmetry::Name to,
  const std::vector<unsigned>& mapping
) {
  /* Creates the list of indices in the target symmetry. Why is this necessary?
   *
   * E.g. An index mapping from linear to T-shaped. The individual
   * symmetry-internal numbering schemes are shown for the symmetry positions.
   *
   *  1  –▶  0
   *  |      |
   * (_)    (_) – 1 (new)
   *  |      |                Linear pos. 0 to
   *  0  –▶  2                pos. 2 in Tshaped
   *                                 |  ┌– Linear pos. 1 to pos. 0 in Tshaped
   *                                 |  |  ┌– This position is new
   * This mapping is represented as {2, 0, 1}.
   *
   * This function writes the indices of original mapping into the target
   * symmetry's indexing scheme.
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
   * Which of these properties are related by target symmetry rotations?
   *
   *
   *     mapping       target indices
   * -----------------------------------
   *    {2, 0, 1}   =>   {1, 2, 0}
   *        ▲                ▲
   *        |                |
   *        X                | C2 rotation in T-shape symmetry
   *        |                |
   *        ▼                ▼
   *    {0, 2, 1}   =>   {0, 2, 1}
   *
   */
  std::vector<unsigned> symmetryPositions (Symmetry::size(to));

  for(unsigned i = 0; i < Symmetry::size(to); ++i) {
    symmetryPositions.at(
      mapping.at(i)
    ) = i;
  }

  return symmetryPositions;
}

DistortionInfo::DistortionInfo(
  std::vector<unsigned> passIndexMapping,
  const double& passAngularDistortion,
  const double& passChiralDistortion
) : indexMapping(std::move(passIndexMapping)),
    angularDistortion(passAngularDistortion),
    chiralDistortion(passChiralDistortion)
{}

std::vector<DistortionInfo> symmetryTransitionMappings(
  const Symmetry::Name symmetryFrom,
  const Symmetry::Name symmetryTo
) {

  /* Symmetries must be adjacent in size (0 = rearrangement,
   * +1 = ligand gain. Ligand loss is a special case where a specific position
   * in the symmetry group is removed, and is not covered here!
   */
  assert(
    (std::set<int> {0, 1}).count(
      static_cast<int>(Symmetry::size(symmetryTo))
      - static_cast<int>(Symmetry::size(symmetryFrom))
    ) == 1
  );

  /* Base idea: We need to go through all possible mappings. In situations where
   * the target symmetry has one more or one fewer ligand, the last index in
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
    Symmetry::size(symmetryFrom),
    Symmetry::size(symmetryTo)
  );

  auto indexMapping = temple::iota<unsigned>(largerSize);

  // Store all distortions calculated
  std::vector<DistortionInfo> distortions;

  // Need to keep track of rotations of mappings to avoid repetition
  std::set<
    std::vector<unsigned>
  > encounteredSymmetryMappings;

  /* Using std::next_permutation generates all possible mappings!
   * do-while is required since the very first mapping must be considered before
   * calling next_permutation for the first time
   */
  do {
    if( // is the mapping new?
      encounteredSymmetryMappings.count(
        applyIndexMapping(
          symmetryTo,
          indexMapping
        )
      ) == 0
    ) {
      /* Add it to the set of possible distortions, calculate angular and
       * chiral distortion involved in the mapping
       */
      distortions.emplace_back(
        indexMapping,
        calculateAngleDistortion(symmetryFrom, symmetryTo, indexMapping),
        calculateChiralDistortion(symmetryFrom, symmetryTo, indexMapping)
      );

      /* Any rotations of the mapping in the target symmetry are equivalent, we
       * do not want to count these as an additional multiplicity, so we
       * generate them and add them to the encountered mappings
       */
      auto allRotations = generateAllRotations(
        symmetryTo,
        applyIndexMapping(
          symmetryTo,
          indexMapping
        )
      );

      encounteredSymmetryMappings.insert(
        allRotations.begin(),
        allRotations.end()
      );
    }
  } while (std::next_permutation(indexMapping.begin(), indexMapping.end()));

  return distortions;
}

std::vector<DistortionInfo> ligandLossTransitionMappings(
  const Symmetry::Name symmetryFrom,
  const Symmetry::Name symmetryTo,
  const unsigned positionInSourceSymmetry
) {
  // Ensure we are dealing with ligand loss
  assert(Symmetry::size(symmetryTo) + 1 == Symmetry::size(symmetryFrom));
  assert(positionInSourceSymmetry < Symmetry::size(symmetryFrom));

  /* Generate the index mapping specific to this position loss in the target
   * symmetry.
   *
   * The possible distortions are determined from the equivalent case of adding
   * the ligand (which was to be deleted) to the smaller symmetry.
   *
   * The deleted index is added to the end. The lowest permutation of the mapping
   * from the smaller to the larger symmetry is then:
   *
   * 1, 2, ..., (pos - 1), (pos + 1), ..., (from_size - 1), pos
   */
  std::vector<unsigned> indexMapping = temple::variadic::concatenate(
    temple::iota<unsigned>(positionInSourceSymmetry),
    detail::range(positionInSourceSymmetry + 1, Symmetry::size(symmetryFrom))
  );

  /* NOTE: From here the algorithm is identical to symmetryTransitionMappings
   * save that symmetryTo and symmetryFrom are swapped in all occasions
   * and that std::next_permutation is only called on the subset excluding the
   * last position (the one that is added / deleted).
   */

  std::vector<DistortionInfo> distortions;

  std::set<
    std::vector<unsigned>
  > encounteredSymmetryMappings;

  do {
    if(encounteredSymmetryMappings.count(indexMapping) == 0) {
      distortions.emplace_back(
        indexMapping,
        calculateAngleDistortion(symmetryTo, symmetryFrom, indexMapping),
        calculateChiralDistortion(symmetryTo, symmetryFrom, indexMapping)
      );

      auto allRotations = generateAllRotations(symmetryTo, indexMapping);

      encounteredSymmetryMappings.insert(
        allRotations.begin(),
        allRotations.end()
      );
    }
  } while (std::next_permutation(indexMapping.begin(), indexMapping.end()));

  return distortions;
}

SymmetryTransitionGroup::SymmetryTransitionGroup(
  std::vector<
    std::vector<unsigned>
  > passIndexMappings,
  const double& passAngleDistortion,
  const double& passChiralDistortion
) : indexMappings(std::move(passIndexMappings)),
    angularDistortion(passAngleDistortion),
    chiralDistortion(passChiralDistortion)
{}

SymmetryTransitionGroup selectBestTransitionMappings(
  const std::vector<DistortionInfo>& distortions
) {
  /* We are interested only in the those transitions that have the very lowest
   * angular distortion, and within that set only the lowest chiral distortion,
   * so we sub-select within the generated set
   */

  double lowestAngularDistortion = temple::min(
    temple::adaptors::transform(
      distortions,
      [](const auto& distortion) -> double {
        return distortion.angularDistortion;
      }
    )
  );

  auto distortionsView = temple::view_filter(
    distortions,
    [&lowestAngularDistortion](const auto& distortion) -> bool {
      return (
        distortion.angularDistortion > (
          lowestAngularDistortion + floatingPointEqualityThreshold
        )
      );
    }
  );

  // And now sub-set further on the lowest chiral distortion
  double lowestChiralDistortion = temple::min(
    temple::adaptors::transform(
      distortionsView,
      [](const auto& distortion) -> double {
        return distortion.chiralDistortion;
      }
    )
  );

  distortionsView.filter(
    [&lowestChiralDistortion](const auto& distortion) -> bool {
      return (
        distortion.chiralDistortion > (
          lowestChiralDistortion + floatingPointEqualityThreshold
        )
      );
    }
  );

  std::vector<
    std::vector<unsigned>
  > mappings;

  for(const auto& distortionInfo : distortionsView) {
    mappings.push_back(distortionInfo.indexMapping);
  }

  return SymmetryTransitionGroup(
    mappings,
    lowestAngularDistortion,
    lowestChiralDistortion
  );
}


unsigned numUnlinkedAssignments(
  const Symmetry::Name symmetry,
  const unsigned nIdenticalLigands
) {
  unsigned count = 1;
  auto indices = detail::range(0u, Symmetry::size(symmetry));

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  std::set<decltype(indices)> rotations;

  auto initialRotations = generateAllRotations(symmetry, indices);

  for(const auto& rotation : initialRotations) {
    rotations.insert(rotation);
  }

  while(std::next_permutation(indices.begin(), indices.end())) {
    if(rotations.count(indices) == 0) {
      auto allRotations = generateAllRotations(symmetry, indices);
      for(const auto& rotation : allRotations) {
        rotations.insert(rotation);
      }

      ++count;
    }
  }

  return count;
}

bool hasMultipleUnlinkedAssignments(
  const Symmetry::Name symmetry,
  const unsigned nIdenticalLigands
) {
  auto indices = detail::range(0u, Symmetry::size(symmetry));

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  std::set<decltype(indices)> rotations;

  auto initialRotations = generateAllRotations(symmetry, indices);

  for(const auto& rotation : initialRotations) {
    rotations.insert(rotation);
  }

  while(std::next_permutation(indices.begin(), indices.end())) {
    if(rotations.count(indices) == 0) {
      return true;
    }
  }

  return false;
}

} // namespace properties

} // namespace Symmetry
