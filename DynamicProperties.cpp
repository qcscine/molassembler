#include "DynamicProperties.h"

#include <Eigen/Dense>

#include "template_magic/MemberFetcher.h"
#include "template_magic/Numeric.h"
#include "template_magic/VectorView.h"

/* TODO
 * - Consider unordered_set for rotation contains, using the hash function from 
 *   constexpr work
 */

namespace Symmetry {

namespace properties {

std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const Symmetry::Name& symmetryName,
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

Eigen::Vector3d getCoordinates(
  const Symmetry::Name& symmetryName,
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
  const Symmetry::Name& from,
  const Symmetry::Name& to,
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
  constexpr unsigned maxDigitsStoreable = ConstexprMagic::Math::floor(
    ConstexprMagic::Math::log10(
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

  for(unsigned i = 0; i < indexList.size(); ++i) {
    if(indexList.at(i) > 9) {
      throw std::logic_error("hashIndexList: Index list contains numbers greater than 9!");
    }
    hash += tenPowers * indexList.at(i);
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
  const Symmetry::Name& from,
  const Symmetry::Name& to,
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
  const Symmetry::Name& symmetryName,
  const std::vector<unsigned>& indices
) {
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
  const Symmetry::Name& to,
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
  const Symmetry::Name& symmetryFrom,
  const Symmetry::Name& symmetryTo
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

  auto indexMapping = detail::iota<unsigned>(largerSize);

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
  const Symmetry::Name& symmetryFrom,
  const Symmetry::Name& symmetryTo,
  const unsigned& positionInSourceSymmetry
) {
  /* TODO this algorithm and it's constexpr counterpart may give the same
   * results, but they are both incorrect. Rotational equivalence of the
   * mapping must be considered in the target symmetry, not in the source
   * symmetry. Merely limiting the permutations considered in the
   * corresponding ligand gain case is insufficient.
   *
   * E.g. square-pyramidal to square-planar. This is considered as
   * square-planar to square-pyramidal ligand gain, where the new ligand is
   * either an equatorial one (0-3 in the sq-py indexing) or apical (4).
   * In the apical case, although there is clearly only one rotationally unique
   * mapping from square-pyramidal to square-planar when the apical ligand is
   * removed (since 1-2-3-4 and 4-3-2-1 are C2' superposable), the inverse case
   * is considered, and there are two possible mappings, in which the
   * positioning of the new apical ligand can be either at the top or bottom of
   * an existing index sequence for square-planar.
   */

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
  std::vector<unsigned> indexMapping = TemplateMagic::concatenate(
    detail::iota<unsigned>(positionInSourceSymmetry),
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

      auto allRotations = generateAllRotations(
        symmetryTo,
        indexMapping
      );

      encounteredSymmetryMappings.insert(
        allRotations.begin(),
        allRotations.end()
      );
    }
  } while (std::next_permutation(indexMapping.begin(), indexMapping.end()));

  return distortions;
}

SymmetryTransitionGroup::SymmetryTransitionGroup(
  std::set<
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

  double lowestAngularDistortion = TemplateMagic::min(
    TemplateMagic::getMember(
      distortions,
      [](const auto& distortion) -> double {
        return distortion.angularDistortion;
      }
    )
  );

  auto distortionsView = TemplateMagic::filter(
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
  double lowestChiralDistortion = TemplateMagic::min(
    TemplateMagic::getMember(
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

  std::set<
    std::vector<unsigned>
  > mappings;

  for(const auto& distortionInfo : distortionsView) {
    mappings.insert(distortionInfo.indexMapping);
  }

  return SymmetryTransitionGroup(
    mappings,
    lowestAngularDistortion,
    lowestChiralDistortion
  );
}


unsigned numUnlinkedAssignments(
  const Symmetry::Name& symmetry,
  const unsigned& nIdenticalLigands
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

} // namespace properties

} // namespace Symmetry
