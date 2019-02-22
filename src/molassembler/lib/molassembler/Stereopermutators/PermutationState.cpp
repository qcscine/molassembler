/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/PermutationState.h"

#include "CyclicPolygons.h"
#include "temple/Adaptors/SequentialPairs.h"
#include "temple/Adaptors/Transform.h"
#include "temple/Functional.h"

#include "molassembler/Cycles.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/Modeling/CommonTrig.h"

namespace Scine {

namespace molassembler {

PermutationState::PermutationState(
  const RankingInformation& ranking,
  const AtomIndex centerAtom,
  const Symmetry::Name symmetry,
  const OuterGraph& graph
) {
  canonicalLigands = canonicalize(ranking.ligandsRanking);
  symbolicCharacters = transferToSymbolicCharacters(canonicalLigands);
  selfReferentialLinks = selfReferentialTransform(
    ranking.links,
    canonicalLigands
  );

  using ModelType = DistanceGeometry::SpatialModel;

  ligandDistances = temple::map(
    ranking.ligands,
    [&](const auto& ligandIndices) -> DistanceGeometry::ValueBounds {
      return ModelType::ligandDistanceFromCenter(
        ligandIndices,
        centerAtom,
        graph
      );
    }
  );

  Cycles etaLessCycles {graph, true};
  coneAngles.reserve(ranking.ligands.size());

  for(unsigned i = 0; i < ranking.ligands.size(); ++i) {
    coneAngles.push_back(
      ModelType::coneAngle(
        ranking.ligands.at(i),
        ligandDistances.at(i),
        graph,
        etaLessCycles
      )
    );
  }

  permutations = stereopermutation::uniquesWithWeights(
    stereopermutation::Stereopermutation {
      symmetry,
      symbolicCharacters,
      selfReferentialLinks
    },
    symmetry,
    false // Do NOT remove trans-spanning linked groups
  );

  // Determine which permutations are feasible and which aren't
  if(
    // Links are present
    !ranking.links.empty()
    // OR there are haptic ligands
    || temple::sum(
      temple::adaptors::transform(
        ranking.ligands,
        [](const auto& ligandIndices) -> unsigned {
          if(ligandIndices.size() == 1) {
            return 0;
          }

          return 1;
        }
      )
    ) > 0
  ) {
    feasiblePermutations.reserve(permutations.stereopermutations.size());
    const unsigned P = permutations.stereopermutations.size();
    for(unsigned i = 0; i < P; ++i) {
      if(
        isNotObviouslyImpossibleStereopermutation(
          permutations.stereopermutations.at(i),
          canonicalLigands,
          coneAngles,
          ranking,
          symmetry,
          graph
        )
      ) {
        feasiblePermutations.push_back(i);
      }
    }
    feasiblePermutations.shrink_to_fit();
  } else {
    feasiblePermutations = temple::iota<unsigned>(permutations.stereopermutations.size());
  }
}

RankingInformation::RankedLigandsType PermutationState::canonicalize(
  RankingInformation::RankedLigandsType rankedLigands
) {
  std::stable_sort(
    rankedLigands.begin(),
    rankedLigands.end(),
    [](const auto& setA, const auto& setB) -> bool {
      // Inverted comparison so that larger sets come first
      return setA.size() > setB.size();
    }
  );

  return rankedLigands;
}

// Transform canonical ranked ligands to canonical characters
std::vector<char> PermutationState::transferToSymbolicCharacters(
  const RankingInformation::RankedLigandsType& canonicalLigands
) {
  std::vector<char> characters;

  char currentChar = 'A';
  for(const auto& equalPrioritySet : canonicalLigands) {
    for(unsigned i = 0; i < equalPrioritySet.size(); ++i) {
      characters.push_back(currentChar);
    }

    ++currentChar;
  }

  return characters;
}

stereopermutation::Stereopermutation::LinksSetType
PermutationState::selfReferentialTransform(
  const std::vector<LinkInformation>& rankingLinks,
  const RankingInformation::RankedLigandsType& canonicalLigands
) {
  stereopermutation::Stereopermutation::LinksSetType links;

  for(const auto& link : rankingLinks) {
    auto getRankedPosition = [&canonicalLigands](const unsigned ligandIndex) -> unsigned {
      unsigned position = 0;
      for(const auto& equalLigandsSet : canonicalLigands) {
        for(const auto& rankedLigandIndex : equalLigandsSet) {
          if(rankedLigandIndex == ligandIndex) {
            return position;
          }

          ++position;
        }
      }

      throw std::logic_error("Ligand index not found in ranked ligands");
    };

    auto a = getRankedPosition(link.indexPair.first),
         b = getRankedPosition(link.indexPair.second);

    links.emplace(
      std::min(a, b),
      std::max(a, b)
    );
  }

  return links;
}

std::vector<unsigned>
PermutationState::generateLigandToSymmetryPositionMap(
  const stereopermutation::Stereopermutation& assignment,
  const RankingInformation::RankedLigandsType& canonicalLigands
) {
  /* We are given an assignment within some symmetry:
   *   Characters: ABADCB
   *   Links: (0, 1), (2, 5)
   * and canonical ligands (ranked sets of ligand indices re-sorted by
   * decreasing set size): {{0, 4}, {2, 1}, {5}, {3}}
   *
   * We need to distribute the ligand indices (from the canonical ligands) to
   * the symmetry positions (defined by the assignment) that they match (via
   * their ranking character).
   *
   * Additionally, we need to ensure that:
   * - AAAAAA {0, 1}, {2, 3}, {4, 5}
   * - AAAAAA {0, 1}, {2, 4}, {3, 5}
   * have different symmetry position maps.
   *
   * Link indices (from assignments) specify symmetry positions that are linked.
   * Since symmetry positions are NOT exchangeable as two ligand indices are
   * that rank equally, we need to distribute linked symmetry positions first,
   * and then distribute the remaining characters afterwards.
   */

  constexpr unsigned placeholder = std::numeric_limits<unsigned>::max();
  const unsigned S = assignment.characters.size();

  std::vector<unsigned> positionMap (S, placeholder);

  /* Process the links */

  /* For every ligand index within the group of ligands of equal priority, we
   * have to keep information on which have been used and which haven't.
   */
  auto usedLists = temple::map(
    canonicalLigands,
    [](const auto& equalPrioritySet) -> std::vector<bool> {
      return std::vector<bool>(equalPrioritySet.size(), false);
    }
  );

  /* Additionally, for each canonical character, a limited set of symmetry
   * positions are available: Those where the passed assignment's characters
   * match the character.
   */
  std::vector<
    std::vector<unsigned>
  > availableSymmetryPositions;

  const char maxChar = temple::max(assignment.characters);
  // For each ranking character
  for(char i = 'A'; i <= maxChar; ++i) {
    std::vector<unsigned> positions;

    // Go through the symmetry positions
    for(unsigned s = 0; s < S; ++s) {
      if(assignment.characters.at(s) == i) {
        positions.push_back(s);
      }
    }

    availableSymmetryPositions.emplace_back(
      std::move(positions)
    );
  }

  // For linked ligands, we need to find a symmetry position to place them
  auto placeAndMark = [&](const unsigned symmetryPosition) {
    char priority = assignment.characters.at(symmetryPosition);

    unsigned countUpToPosition = std::count(
      std::begin(assignment.characters),
      std::begin(assignment.characters) + symmetryPosition,
      priority
    );

    unsigned correspondingLigand = canonicalLigands.at(priority - 'A').at(countUpToPosition);

    if(positionMap.at(correspondingLigand) == placeholder) {
      unsigned newSymmetryPosition = availableSymmetryPositions.at(priority - 'A').front();

      positionMap.at(correspondingLigand) = newSymmetryPosition;

      availableSymmetryPositions.at(priority - 'A').erase(
        availableSymmetryPositions.at(priority - 'A').begin()
      );

      usedLists.at(priority - 'A').at(countUpToPosition) = true;
    }
  };

  // Place all linked indices
  for(const auto& link: assignment.links) {
    placeAndMark(link.first);
    placeAndMark(link.second);
  }

  // Next, process all characters
  for(const auto& priorityChar : assignment.characters) {
    // Get an unused ligand index for that priority
    const auto unusedIndexIter = std::find(
      std::begin(usedLists.at(priorityChar - 'A')),
      std::end(usedLists.at(priorityChar - 'A')),
      false
    );

    if(unusedIndexIter != std::end(usedLists.at(priorityChar - 'A'))) {
      const unsigned correspondingLigand = canonicalLigands.at(priorityChar - 'A').at(
        unusedIndexIter - std::begin(usedLists.at(priorityChar - 'A'))
      );

      assert(positionMap.at(correspondingLigand) == placeholder);

      const unsigned symmetryPosition = availableSymmetryPositions.at(priorityChar - 'A').front();

      availableSymmetryPositions.at(priorityChar - 'A').erase(
        std::begin(availableSymmetryPositions.at(priorityChar - 'A'))
      );

      positionMap.at(correspondingLigand) = symmetryPosition;

      *unusedIndexIter = true;
    }
  }

  // Ensure no symmetry positions are marked with placeholders
  assert(
    temple::all_of(
      positionMap,
      [](const unsigned symmetryPosition) -> bool {
        return symmetryPosition != placeholder;
      }
    ) && "A symmetry position is still marked with a placeholder!"
  );

  return positionMap;
}

std::vector<unsigned> PermutationState::generateSymmetryPositionToLigandMap(
  const stereopermutation::Stereopermutation& assignment,
  const RankingInformation::RankedLigandsType& canonicalLigands
) {
  auto base = generateLigandToSymmetryPositionMap(assignment, canonicalLigands);

  std::vector<unsigned> inverseMap (base.size());

  for(unsigned i = 0; i < base.size(); ++i) {
    inverseMap.at(base.at(i)) = i;
  }

  return inverseMap;
}

std::vector<char> PermutationState::makeStereopermutationCharacters(
  const RankingInformation::RankedLigandsType& canonicalLigands,
  const std::vector<char>& canonicalStereopermutationCharacters,
  const std::vector<unsigned>& ligandsAtSymmetryPositions
) {
  // Replace the ligand indices by their new ranking characters
  std::vector<unsigned> flattenedIndices;
  for(const auto& equalPrioritySet : canonicalLigands) {
    for(const auto& index : equalPrioritySet) {
      flattenedIndices.push_back(index);
    }
  }

  std::vector<char> newStereopermutationCharacters;

  for(const auto& index : ligandsAtSymmetryPositions) {
    const auto findIter = std::find(
      flattenedIndices.begin(),
      flattenedIndices.end(),
      index
    );

    newStereopermutationCharacters.push_back(
      canonicalStereopermutationCharacters.at(
        findIter - flattenedIndices.begin()
      )
    );
  }

  return newStereopermutationCharacters;
}

/* WARNING: This has to be a copy-initialized optional. Don't change it, unless
 * you want to sift through -fsanitize=address output to find the bug in the
 * optional propagation. It's not safe to return a reference to within a
 * temporary object where this is used.
 */
boost::optional<std::vector<unsigned>> PermutationState::getIndexMapping(
  const Symmetry::properties::SymmetryTransitionGroup& mappingsGroup,
  const ChiralStatePreservation& preservationOption
) {
  if(mappingsGroup.indexMappings.empty()) {
    return boost::none;
  }

  if(
    preservationOption == ChiralStatePreservation::EffortlessAndUnique
    && mappingsGroup.indexMappings.size() == 1
    && mappingsGroup.angularDistortion <= 0.2
  ) {
    return mappingsGroup.indexMappings.front();
  }

  if(
    preservationOption == ChiralStatePreservation::Unique
    && mappingsGroup.indexMappings.size() == 1
  ) {
    return mappingsGroup.indexMappings.front();
  }

  if(preservationOption == ChiralStatePreservation::RandomFromMultipleBest) {
    return mappingsGroup.indexMappings.at(
      temple::random::getSingle<unsigned>(
        0,
        mappingsGroup.indexMappings.size() - 1,
        randomnessEngine()
      )
    );
  }

  return boost::none;
}

bool PermutationState::isNotObviouslyImpossibleStereopermutation(
  const stereopermutation::Stereopermutation& assignment,
  const RankingInformation::RankedLigandsType& canonicalLigands,
  const ConeAngleType& coneAngles,
  const RankingInformation& ranking,
  const Symmetry::Name symmetry,
  const OuterGraph& graph
) {
  auto symmetryPositionMap = generateLigandToSymmetryPositionMap(
    assignment,
    canonicalLigands
  );

  // Check if any haptic ligand cones intersect
  const unsigned L = ranking.ligands.size();
  for(unsigned ligandI = 0; ligandI < L - 1; ++ligandI) {
    for(unsigned ligandJ = ligandI + 1; ligandJ < L; ++ligandJ) {
      // Do not test cone angles if no angle could be calculated
      if(!coneAngles.at(ligandI) || !coneAngles.at(ligandJ)) {
        continue;
      }

      double symmetryAngle = Symmetry::angleFunction(symmetry)(
        symmetryPositionMap.at(ligandI),
        symmetryPositionMap.at(ligandJ)
      );

      /* A haptic steropermutation of ligands is only obviously impossible if
       * the haptic ligands have no spatial freedom to arrange in a fashion
       * that does not overlap.
       */
      if(
        (
          symmetryAngle
          - coneAngles.at(ligandI).value().lower
          - coneAngles.at(ligandJ).value().lower
        ) < 0
      ) {
        return false;
      }
    }
  }

  /* Idea: An assignment is unfeasible if any link's cycle cannot be realized
   * as a flat cyclic polygon, in which the edges from the central atom are
   * merged using the joint angle calculable from the assignment and symmetry.
   *
   * The algorithm below is explained in detail in
   * documents/denticity_feasibility/.
   */
  for(const auto& link : ranking.links) {
    // Ignore cycles of size 3
    if(link.cycleSequence.size() == 4) {
      continue;
    }

    // Perform no checks if, for either of the ligands, no cone angle could be calculated
    if(!coneAngles.at(link.indexPair.first) || !coneAngles.at(link.indexPair.second)) {
      continue;
    }

    const DistanceGeometry::ValueBounds ligandIConeAngle = coneAngles.at(link.indexPair.first).value();
    const DistanceGeometry::ValueBounds ligandJConeAngle = coneAngles.at(link.indexPair.second).value();

    const double symmetryAngle = Symmetry::angleFunction(symmetry)(
      symmetryPositionMap.at(link.indexPair.first),
      symmetryPositionMap.at(link.indexPair.second)
    );

    /* A link across haptic ligands is only obviously impossible if it is
     * impossible in the best case scenario. In this case, especially for alpha,
     * ligand bridge links must be possible only in the best case spatial
     * arrangement for the haptic ligand link to be possible. That means
     * subtracting the upper bound of the respective cone angles.
     */
    const double alpha = std::max(
      0.0,
      symmetryAngle - ligandIConeAngle.upper - ligandJConeAngle.upper
    );

    auto edgeLengths = temple::map(
      temple::adaptors::sequentialPairs(link.cycleSequence),
      [&](const auto& i, const auto& j) -> double {
        return Bond::calculateBondDistance(
          graph.elementType(i),
          graph.elementType(j),
          graph.bondType(BondIndex {i, j})
        );
      }
    );

    const double a = edgeLengths.front();
    const double b = edgeLengths.back();
    const double c = CommonTrig::lawOfCosines(a, b, alpha);

    edgeLengths.front() = c;
    edgeLengths.erase(edgeLengths.end() - 1);

    // Quick escape: If the cyclic polygon isn't even constructible, fail
    if(!CyclicPolygons::exists(edgeLengths)) {
      return false;
    }

    /* Test that no atom in cyclic polygon except binding sites in binding
     * distance to central atom.
     */

    const auto phis = CyclicPolygons::internalAngles(edgeLengths);

    const double d1 = CommonTrig::lawOfCosines(
      a,
      edgeLengths.at(1), // i
      phis.at(0) + alpha
    );

    if(
      d1 <= Bond::calculateBondDistance(
        graph.elementType(link.cycleSequence.at(0)),
        // 0 is the central index, 1 is the first ligand, etc.
        graph.elementType(link.cycleSequence.at(2)),
        BondType::Single
      )
    ) {
      return false;
    }

    std::vector<double> distances {a, d1};
    std::vector<double> deltas {};

    for(unsigned i = 1; i < phis.size() - 2; ++i) {
      deltas.push_back(
        CommonTrig::lawOfCosinesAngle(
          edgeLengths.at(i),
          *(distances.end() - 1),
          *(distances.end() - 2)
        )
      );

      distances.push_back(
        CommonTrig::lawOfCosines(
          distances.back(),
          edgeLengths.at(i + 1),
          phis.at(i) - deltas.back()
        )
      );

      if(
        distances.back() <= Bond::calculateBondDistance(
          graph.elementType(link.cycleSequence.at(0)),
          // 0 is the central index, 1 is the first ligand, etc.
          graph.elementType(link.cycleSequence.at(i + 2)),
          BondType::Single
        )
      ) {
        return false;
      }
    }
  }

  return true;
}

} // namespace molassembler

} // namespace Scine
