#include "temple/constexpr/ConsecutiveCompare.h"
#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Optionals.h"
#include "temple/Random.h"

#include "chemical_symmetries/Properties.h"
#include "CyclicPolygons.h"

#include "boost/range/combine.hpp"
#include "boost/range/algorithm/remove_if.hpp"

#include "Molecule.h"
#include "BuildTypeSwitch.h"
#include "CNStereocenter.h"
#include "CommonTrig.h"
#include "DelibHelpers.h"
#include "Log.h"
#include "StdlibTypeAlgorithms.h"

#include <iomanip>

namespace molassembler {

namespace Stereocenters {

namespace adhesive {

RankingInformation::RankedType ligandRanking(
  const RankingInformation::RankedType& sortedSubstituents,
  const RankingInformation::RankedType& ligands
) {
  // TODO templify the find() call
  auto ligandCoefficients = temple::map(
    ligands,
    [&sortedSubstituents](const auto& ligandSet) -> unsigned {
      return temple::sum(
        temple::map(
          ligandSet,
          [&sortedSubstituents](const auto& ligandIndex) -> unsigned {
            return temple::find_if(
              sortedSubstituents,
              [&ligandIndex](const auto& equalPrioritySet) -> bool {
                return temple::find(equalPrioritySet, ligandIndex) != equalPrioritySet.end();
              }
            ) - sortedSubstituents.begin();
          }
        )
      );
    }
  );

  std::vector<unsigned> sortPermutation (ligands.size());
  std::iota(sortPermutation.begin(), sortPermutation.end(), 0);
  std::sort(
    sortPermutation.begin(),
    sortPermutation.end(),
    [&ligandCoefficients](const auto& a, const auto& b) -> bool {
      return ligandCoefficients.at(a) < ligandCoefficients.at(b);
    }
  );

  /* Since the sort does not group ligands of equal coefficients, we have to do
   * that in a separate step:
   */

  RankingInformation::RankedType rankedLigands;
  for(const auto& ligandPosition : sortPermutation) {
    if(rankedLigands.size() == 0) {
      rankedLigands.emplace_back();
      rankedLigands.back().push_back(ligandPosition);
    } else {
      if(
        ligandCoefficients.at(ligandPosition)
        == ligandCoefficients.at(rankedLigands.back().back())
      ) {
        rankedLigands.back().push_back(ligandPosition);
      } else {
        rankedLigands.emplace_back();
        rankedLigands.back().push_back(ligandPosition);
      }
    }
  }

  return rankedLigands;
}

/* Can use output of ligandRanking directly to get canonical characters.
 * NOTE: this is the same as in glue
 */
std::vector<char> canonicalCharacters(
  const RankingInformation::RankedType& rankedLigands
) {
  std::vector<char> characters;

  char currentChar = 'A';
  for(const auto& equalPrioritySet : rankedLigands) {
    for(unsigned i = 0; i < equalPrioritySet.size(); ++i) {
      characters.push_back(currentChar);
    }

    ++currentChar;
  }

  return characters;
}

stereopermutation::Stereopermutation::LinksSetType canonicalLinks(
  const RankingInformation::RankedType& ligands,
  const RankingInformation::RankedType& rankedLigands,
  const RankingInformation::LinksType& rankingLinks
) {
  stereopermutation::Stereopermutation::LinksSetType links;

  for(const auto& link : rankingLinks) {
    auto getLigandIndex = [&ligands](const AtomIndexType& soughtIndex) -> AtomIndexType {
      return std::find_if(
        ligands.begin(),
        ligands.end(),
        [&soughtIndex](const auto& ligandMembers) -> bool {
          return std::find(
            ligandMembers.begin(),
            ligandMembers.end(),
            soughtIndex
          ) != ligandMembers.end();
        }
      ) - ligands.begin();
    };

    auto getRankedPosition = [&rankedLigands](const AtomIndexType& ligandIndex) -> unsigned {
      unsigned position = 0;
      for(const auto& equalLigandsSet : rankedLigands) {
        for(const auto& rankedLigandIndex : equalLigandsSet) {
          if(rankedLigandIndex == ligandIndex) {
            return position;
          }

          ++position;
        }
      }

      throw std::logic_error("Ligand index not found in ranked ligands");
    };

    auto a = getRankedPosition(getLigandIndex(link.indexPair.first)),
         b = getRankedPosition(getLigandIndex(link.indexPair.second));

    links.emplace(
      std::min(a, b),
      std::max(a, b)
    );
  }

  return links;
}

} // namespace adhesive

namespace glue {

RankingInformation::RankedType canonicalize(
  RankingInformation::RankedType sortedSubstituents
) {
  /* Stable sort so that ranking information is preserved while moving around
   * sets of different sizes
   */
  std::stable_sort(
    sortedSubstituents.begin(),
    sortedSubstituents.end(),
    [](const auto& setA, const auto& setB) -> bool {
      // Inverted comparison so that larger sets come first
      return setA.size() > setB.size();
    }
  );

  return sortedSubstituents;
}

std::vector<char> makeCanonicalCharacters(
  const RankingInformation::RankedType& canonicalizedSubstituents
) {
  std::vector<char> characters;

  char currentChar = 'A';
  for(const auto& equalPrioritySet : canonicalizedSubstituents) {
    for(unsigned i = 0; i < equalPrioritySet.size(); ++i) {
      characters.push_back(currentChar);
    }

    ++currentChar;
  }

  return characters;
}

stereopermutation::Stereopermutation::LinksSetType makeLinksSelfReferential(
  const RankingInformation::RankedType& canonicalizedSubstituents,
  const RankingInformation::LinksType& rankingLinks
) {
  // TODO there's definitely a better way to do this

  // Flatten the sorted list of indices
  std::vector<unsigned> sortedIndices;
  for(const auto& equalPrioritySet : canonicalizedSubstituents) {
    for(const auto& index : equalPrioritySet) {
      sortedIndices.push_back(index);
    }
  }

  /* Have a sequence of indices ranked by priority low to high
   * and a set of pairs with the same indices, but we want a set that is
   * self-referential, i.e. refers to the indices in sortedIndices, e.g.:
   *
   * sortedIndices: vec {7, 8, 3, 1}
   * linkedPairsSet: set { pair {7, 8}, pair {7, 3} }
   *
   * -> links: set { pair{0, 1}, pair{0, 2} }
   */

  stereopermutation::Stereopermutation::LinksSetType links;

  auto findIndexInSorted = [&](const AtomIndexType& index) -> unsigned {
    auto findIter = std::find(
      sortedIndices.begin(),
      sortedIndices.end(),
      index
    );

    assert(findIter != sortedIndices.end());

    return findIter - sortedIndices.begin();
  };

  for(const auto& linkPair : rankingLinks) {
    auto a = findIndexInSorted(linkPair.indexPair.first);
    auto b = findIndexInSorted(linkPair.indexPair.second);

    links.emplace(
      std::min(a, b),
      std::max(a, b)
    );
  }

  return links;
}

std::map<AtomIndexType, unsigned> makeSymmetryPositionMap(
  const stereopermutation::Stereopermutation& assignment,
  const RankingInformation::RankedType& canonicalizedSubstituents
) {
  /* TODO NO Stereopermutation link information is used here!
   *
   * Need to ensure that
   * AAAAAA {0, 1}, {2, 3}, {4, 5}
   * and AAAAAA {0, 1}, {2, 4}, {3, 5}
   * are different!
   *
   * AABCAD {0,1}, {0, 4}, {1, 4} has to work too
   * - cannot use index of symmetry position to fetch atom index within equal priority groups
   * - must be two-stage algorithm that fixes linked symmetry positions first, then distributes
   *   the remaining
   */

  std::map<AtomIndexType, unsigned> positionMap;

  /* First, process the links.
   *
   * For every atom index within the group of indices of equal priority, we
   * have to keep information on which have been used and which haven't.
   */
  auto usedLists = temple::map(
    canonicalizedSubstituents,
    [](const auto& equalPrioritySet) -> std::vector<bool> {
      return std::vector<bool>(equalPrioritySet.size(), false);
    }
  );

  std::vector<
    std::vector<unsigned>
  > availableSymmetryPositions;

  for(char i = 'A'; i <= temple::max(assignment.characters); ++i) {
    std::vector<unsigned> positions;

    for(unsigned s = 0; s < assignment.characters.size(); ++s) {
      if(assignment.characters.at(s) == i) {
        positions.push_back(s);
      }
    }

    availableSymmetryPositions.emplace_back(
      std::move(positions)
    );
  }

  auto placeAndMark = [&](const unsigned& symmetryPositionIndex) {
    char priority = assignment.characters.at(symmetryPositionIndex);

    unsigned countUpToPosition = std::count(
      assignment.characters.begin(),
      assignment.characters.begin() + symmetryPositionIndex,
      priority
    );

    AtomIndexType correspondingAtom = canonicalizedSubstituents.at(priority - 'A').at(countUpToPosition);

    if(positionMap.count(correspondingAtom) == 0) {
      unsigned newSymmetryPosition = availableSymmetryPositions.at(priority - 'A').front();

      positionMap.emplace(
        correspondingAtom,
        newSymmetryPosition
      );

      availableSymmetryPositions.at(priority - 'A').erase(
        availableSymmetryPositions.at(priority - 'A').begin()
      );

      usedLists.at(priority - 'A').at(countUpToPosition) = true;
    }
  };

  for(const auto& link: assignment.links) {
    placeAndMark(link.first);
    placeAndMark(link.second);
  }

  // Next, process the characters
  for(const auto& priorityChar : assignment.characters) {
    // Get an unused atom index for that priority
    auto unusedIndexIter = std::find(
      usedLists.at(priorityChar - 'A').begin(),
      usedLists.at(priorityChar - 'A').end(),
      false
    );

    if(unusedIndexIter != usedLists.at(priorityChar - 'A').end()) {
      AtomIndexType correspondingAtom = canonicalizedSubstituents.at(priorityChar - 'A').at(
        unusedIndexIter - usedLists.at(priorityChar - 'A').begin()
      );

      assert(positionMap.count(correspondingAtom) == 0);

      unsigned symmetryPosition = availableSymmetryPositions.at(priorityChar - 'A').front();

      availableSymmetryPositions.at(priorityChar - 'A').erase(
        availableSymmetryPositions.at(priorityChar - 'A').begin()
      );

      positionMap.emplace(
        correspondingAtom,
        symmetryPosition
      );

      *unusedIndexIter = true;
    }
  }

  return positionMap;
}

// Just returns a flat, inverted map from makeSymmetryPositionMap
// TODO remove?
std::vector<AtomIndexType> mapToSymmetryPositions(
  const stereopermutation::Stereopermutation& assignment,
  const RankingInformation::RankedType& canonicalizedSubstituents
) {
  auto base = makeSymmetryPositionMap(assignment, canonicalizedSubstituents);
  std::vector<AtomIndexType> flatMap;
  flatMap.resize(base.size());

  for(const auto& iterPair : base) {
    flatMap.at(iterPair.second) = iterPair.first;
  }

  return flatMap;
}

std::vector<char> makeStereopermutationCharacters(
  const RankingInformation::RankedType& canonicalizedSubstituents,
  const std::vector<char>& canonicalizedStereopermutationCharacters,
  const std::vector<AtomIndexType>& atomsAtSymmetryPositions
) {
  // Replace the atom indices by their new ranking characters
  std::vector<unsigned> flattenedIndices;
  for(const auto& equalPrioritySet : canonicalizedSubstituents) {
    for(const auto& index : equalPrioritySet) {
      flattenedIndices.push_back(index);
    }
  }

  std::vector<char> newStereopermutationCharacters;

  for(const auto& index : atomsAtSymmetryPositions) {
    const auto findIter = std::find(
      flattenedIndices.begin(),
      flattenedIndices.end(),
      index
    );

    newStereopermutationCharacters.push_back(
      canonicalizedStereopermutationCharacters.at(
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
boost::optional<std::vector<unsigned>> getIndexMapping(
  const Symmetry::properties::SymmetryTransitionGroup& mappingsGroup,
  const ChiralStatePreservation& preservationOption
) {
  if(mappingsGroup.indexMappings.size() == 0) {
    return boost::none;
  }

  if(preservationOption == ChiralStatePreservation::EffortlessAndUnique) {
    if(
      mappingsGroup.indexMappings.size() == 1
      && mappingsGroup.angularDistortion <= 0.2
    ) {
      return mappingsGroup.indexMappings.front();
    }

    return boost::none;
  }

  if(preservationOption == ChiralStatePreservation::Unique) {
    if(mappingsGroup.indexMappings.size() == 1) {
      return mappingsGroup.indexMappings.front();
    }

    return boost::none;
  }

  if(preservationOption == ChiralStatePreservation::RandomFromMultipleBest) {
    return mappingsGroup.indexMappings.at(
      temple::random.getSingle<unsigned>(
        0,
        mappingsGroup.indexMappings.size() - 1
      )
    );
  }

  return boost::none;
}

} // namespace glue

void CNStereocenter::_removeImpossibleStereopermutations(
  stereopermutation::StereopermutationsWithWeights& data,
  const Molecule& molecule
) {
  auto toRemove = temple::map(
    data.assignments,
    [&](const auto& assignment) -> bool {
      return !isFeasibleStereopermutation(assignment, molecule, _symmetry, _ranking);
    }
  );

  data.assignments.erase(
    std::remove_if(
      data.assignments.begin(),
      data.assignments.end(),
      [&](const auto& assignment) -> bool {
        // For this wonderful pointer arithmetic, see https://stackoverflow.com/a/23123481
        return toRemove.at(&assignment - &*data.assignments.begin());
      }
    ),
    data.assignments.end()
  );

  data.weights.erase(
    std::remove_if(
      data.weights.begin(),
      data.weights.end(),
      [&](const auto& weight) -> bool {
        return toRemove.at(&weight - &*data.weights.begin());
      }
    ),
    data.weights.end()
  );
}


/* Static functions */
bool CNStereocenter::isFeasibleStereopermutation(
  const StereopermutationType& assignment,
  const Molecule& molecule,
  const Symmetry::Name& symmetry,
  const RankingInformation& ranking
) {
  if(ranking.links.size() == 0) {
    return true;
  }

  /* Idea: An assignment is unfeasible if any link's cycle cannot be realized
   * as a flat cyclic polygon, in which the edges from the central atom are
   * merged using the joint angle calculable from the assignment and symmetry.
   *
   * The algorithm below is explained in detail in
   * documents/denticity_feasibility/.
   */
  auto symmetryPositionMap = glue::makeSymmetryPositionMap(
    assignment,
    glue::canonicalize(ranking.sortedSubstituents)
  );

  for(const auto& link : ranking.links) {
    // Ignore cycles of size 3
    if(link.cycleSequence.size() == 4) {
      continue;
    }

    double sigma = Symmetry::angleFunction(symmetry)(
      symmetryPositionMap.at(link.indexPair.first),
      symmetryPositionMap.at(link.indexPair.second)
    );

    auto edgeLengths = temple::mapSequentialPairs(
      link.cycleSequence,
      [&](const auto& i, const auto& j) -> double {
        return Bond::calculateBondDistance(
          molecule.getElementType(i),
          molecule.getElementType(j),
          molecule.getBondType(i, j).value()
        );
      }
    );

    double a = edgeLengths.front();
    double b = edgeLengths.back();
    double c = CommonTrig::lawOfCosines(a, b, sigma);

    edgeLengths.front() = c;
    edgeLengths.erase(edgeLengths.end() - 1);

    // Quick escape: If the cyclic polygon isn't even constructible, fail
    if(!CyclicPolygons::exists(edgeLengths)) {
      return false;
    }

    /* Test that no atom in cyclic polygon except binding sites in binding
     * distance to central atom.
     */

    auto phis = CyclicPolygons::internalAngles(edgeLengths);

    double alpha = CommonTrig::lawOfCosinesAngle(a, c, b);

    double d1 = CommonTrig::lawOfCosines(
      a,
      edgeLengths.at(1), // i
      phis.at(0) + alpha
    );

    if(
      d1 <= Bond::calculateBondDistance(
        molecule.getElementType(link.cycleSequence.at(0)),
        // 0 is the central index, 1 is the first ligand
        molecule.getElementType(link.cycleSequence.at(2)),
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
          molecule.getElementType(link.cycleSequence.at(0)),
          // 0 is the central index, 1 is the first ligand
          molecule.getElementType(link.cycleSequence.at(i + 2)),
          BondType::Single
        )
      ) {
        return false;
      }
    }
  }

  return true;
}

/* Constructors */
CNStereocenter::CNStereocenter(
  const Molecule& molecule,
  // The symmetry of this Stereocenter
  const Symmetry::Name& symmetry,
  // The atom this Stereocenter is centered on
  const AtomIndexType& centerAtom,
  // Ranking information of substituents
  const RankingInformation& ranking
) : _ranking {ranking},
    _centerAtom {centerAtom},
    _symmetry {symmetry},
    _assignmentOption {boost::none}
{
  // Canonicalize the ranking's substituents
  _ranking.sortedSubstituents = glue::canonicalize(_ranking.sortedSubstituents);

  // Generate the set of unique assignments possible here
  _uniqueStereopermutationsCache = stereopermutation::uniqueStereopermutationsWithWeights(
    StereopermutationType(
      symmetry,
      glue::makeCanonicalCharacters(_ranking.sortedSubstituents),
      glue::makeLinksSelfReferential(_ranking.sortedSubstituents, _ranking.links)
    ),
    symmetry,
    false // Do NOT remove trans-spanning linked groups
  );

  _removeImpossibleStereopermutations(_uniqueStereopermutationsCache, molecule);
}

/* Modification */
void CNStereocenter::addSubstituent(
  const Molecule& molecule,
  const AtomIndexType& newSubstituentIndex,
  const RankingInformation& newRanking,
  const Symmetry::Name& newSymmetry,
  const ChiralStatePreservation& preservationOption
) {
  auto canonicalizedSubstituents = glue::canonicalize(newRanking.sortedSubstituents);
  auto canonicalizedCharacters = glue::makeCanonicalCharacters(canonicalizedSubstituents);
  auto newLinks = glue::makeLinksSelfReferential(canonicalizedSubstituents, newRanking.links);
  auto newStereopermutations = stereopermutation::uniqueStereopermutationsWithWeights(
    StereopermutationType {
      newSymmetry,
      canonicalizedCharacters,
      newLinks
    },
    newSymmetry
  );

  _removeImpossibleStereopermutations(newStereopermutations, molecule);

  boost::optional<unsigned> newStereopermutation = boost::none;

  /* Does the current stereocenter carry chiral information?
   * If so, try to get a mapping to the new symmetry
   * If that returns a Some, try to get a mapping by preservationOption policy
   *
   * If any of these steps returns boost::none, the whole expression is
   * boost::none.
   */
  auto suitableMappingOption = _assignmentOption
    | temple::callIfSome(Symmetry::getMapping, _symmetry, newSymmetry, boost::none)
    | temple::callIfSome(
        glue::getIndexMapping,
        temple::ANS, // Inserts getMapping's optional value here
        preservationOption
      );

  if(suitableMappingOption) {
    /* So now we must transfer the current assignment into the new symmetry
     * and search for it in the set of uniques.
     */
    const auto& symmetryMapping = suitableMappingOption.value();

    // Transform current assignment from characters to indices
    const auto& currentStereopermutation = _uniqueStereopermutationsCache.assignments.at(
      _assignmentOption.value()
    );

    std::vector<AtomIndexType> atomsAtSmallerSymmetryPositions = glue::mapToSymmetryPositions(
      currentStereopermutation,
      _ranking.sortedSubstituents
    );

    atomsAtSmallerSymmetryPositions.push_back(newSubstituentIndex);

    // Transfer indices from smaller symmetry to larger
    std::vector<AtomIndexType> atomsAtLargerSymmetryPositions (Symmetry::size(newSymmetry));
    for(unsigned i = 0; i < Symmetry::size(newSymmetry); ++i) {
      atomsAtLargerSymmetryPositions.at(
        symmetryMapping.at(i)
      ) = atomsAtSmallerSymmetryPositions.at(i);
    }

    // Get character representation in new symmetry
    std::vector<char> charactersInLargerSymmetry = glue::makeStereopermutationCharacters(
      canonicalizedSubstituents,
      canonicalizedCharacters,
      atomsAtLargerSymmetryPositions
    );

    // Construct an assignment from it
    auto trialStereopermutation = stereopermutation::Stereopermutation(
      newSymmetry,
      charactersInLargerSymmetry,
      newLinks
    );

    // Generate the rotational equivalents
    auto allTrialRotations = trialStereopermutation.generateAllRotations(newSymmetry);

    // Search for a match from the vector of uniques
    for(unsigned i = 0; i < newStereopermutations.assignments.size(); ++i) {
      if(allTrialRotations.count(newStereopermutations.assignments.at(i)) > 0) {
        newStereopermutation = i;
        break;
      }
    }
  }

  /* Since either there is no steric information present or no unique
   * minimal-effort mapping exists to the new symmetry, it is impossible to
   * unambiguously choose a new assignment. In those cases, newStereopermutation
   * remains boost::none, and this stereocenter loses any chiral information it
   * may have had.
   */

  // Overwrite class state
  _ranking = newRanking;
  _ranking.sortedSubstituents = canonicalizedSubstituents;
  _symmetry = newSymmetry;
  _uniqueStereopermutationsCache = newStereopermutations;
  assign(newStereopermutation);
}

void CNStereocenter::assign(const boost::optional<unsigned>& assignment) {
  if(assignment) {
    assert(assignment.value() < _uniqueStereopermutationsCache.assignments.size());
  }

  // Store current assignment
  _assignmentOption = assignment;

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndexType -> unsigned).
   */
  if(assignment) {
    _symmetryPositionMapCache = glue::makeSymmetryPositionMap(
      _uniqueStereopermutationsCache.assignments.at(assignment.value()),
      _ranking.sortedSubstituents
    );
  } else { // Wipe the map
    _symmetryPositionMapCache.clear();
  }
}

void CNStereocenter::assignRandom() {
  assign(
    temple::random.pickDiscrete(_uniqueStereopermutationsCache.weights)
  );
}

void CNStereocenter::propagateGraphChange(
  const Molecule& molecule,
  const RankingInformation& newRanking
) {
  auto newCanonicalizedSubstituents = glue::canonicalize(newRanking.sortedSubstituents);

  // Has anything about the new ranking changed relative to the old?
  if(
    newCanonicalizedSubstituents != _ranking.sortedSubstituents
    || newRanking.links != _ranking.links
  ) {
    auto newCharacters = glue::makeCanonicalCharacters(newCanonicalizedSubstituents);
    auto newLinks = glue::makeLinksSelfReferential(
      newCanonicalizedSubstituents,
      newRanking.links
    );
    auto newStereopermutations = stereopermutation::uniqueStereopermutationsWithWeights(
      StereopermutationType(_symmetry, newCharacters, newLinks),
      _symmetry,
      false // do NOT remove trans-spanning ligand groups
    );

    _removeImpossibleStereopermutations(newStereopermutations, molecule);

    boost::optional<unsigned> newStereopermutation = boost::none;

    /* Before we overwrite class state, we need to figure out which assignment
     * in the new set of assignments corresponds to the one we have now.
     * This is only necessary in the case that the stereocenter is currently
     * assigned and only possible if the new number of assignments is smaller or
     * equal to the amount we have currently.
     */
    if(
      _assignmentOption
      && newStereopermutations.assignments.size() <= _uniqueStereopermutationsCache.assignments.size()
    ) {
      const auto& currentStereopermutation = _uniqueStereopermutationsCache.assignments.at(
        _assignmentOption.value()
      );

      // Replace the characters by their corresponding indices from the old ranking
      std::vector<AtomIndexType> atomsAtSymmetryPositions = glue::mapToSymmetryPositions(
        currentStereopermutation,
        _ranking.sortedSubstituents
      );

      // Replace the atom indices by their new ranking characters
      std::vector<char> newStereopermutationCharacters = glue::makeStereopermutationCharacters(
        newCanonicalizedSubstituents,
        newCharacters,
        atomsAtSymmetryPositions
      );

      // Create a new assignment with those characters
      auto trialStereopermutation = stereopermutation::Stereopermutation(
        _symmetry,
        newStereopermutationCharacters,
        newLinks
      );

      // Generate all rotations of this trial assignment
      auto allTrialRotations = trialStereopermutation.generateAllRotations(_symmetry);

      // Find out which of the new assignments has a rotational equivalent
      for(unsigned i = 0; i < newStereopermutations.assignments.size(); ++i) {
        if(allTrialRotations.count(newStereopermutations.assignments.at(i)) > 0) {
          newStereopermutation = i;
          break;
        }
      }
    }

    // Overwrite the class state
    _ranking = newRanking;
    _ranking.sortedSubstituents = newCanonicalizedSubstituents;
    _uniqueStereopermutationsCache = newStereopermutations;
    assign(newStereopermutation);
  }
}

void CNStereocenter::propagateVertexRemoval(const AtomIndexType& removedIndex) {
  auto updateIndexInplace = [&removedIndex](AtomIndexType& index) -> void {
    if(index > removedIndex) {
      --index;
    } else if(index == removedIndex) {
      index = std::numeric_limits<AtomIndexType>::max();
    }
  };

  auto updateIndex = [&removedIndex](const AtomIndexType& index) -> AtomIndexType {
    if(index > removedIndex) {
      return index - 1;
    }

    if(index == removedIndex) {
      return std::numeric_limits<AtomIndexType>::max();
    }

    return index;
  };


  for(auto& equalPrioritySet : _ranking.sortedSubstituents) {
    for(auto& index : equalPrioritySet) {
      updateIndexInplace(index);
    }
  }

  for(auto& link : _ranking.links) {
    link.indexPair = {
      std::min(updateIndex(link.indexPair.first), updateIndex(link.indexPair.second)),
      std::max(updateIndex(link.indexPair.first), updateIndex(link.indexPair.second))
    };

    link.cycleSequence = temple::map(
      link.cycleSequence,
      updateIndex
    );
  }

  updateIndexInplace(_centerAtom);

  if(_assignmentOption) {
    std::map<AtomIndexType, unsigned> newSymmetryPositionMap;

    for(const auto& iterPair : _symmetryPositionMapCache) {
      const auto& atomIndex = iterPair.first;
      const auto& symmetryPosition = iterPair.second;

      newSymmetryPositionMap.emplace(
        updateIndex(atomIndex),
        symmetryPosition
      );
    }

    _symmetryPositionMapCache = std::move(newSymmetryPositionMap);
  }
}

void CNStereocenter::removeSubstituent(
  const Molecule& molecule,
  const AtomIndexType& which,
  const RankingInformation& newRanking,
  const Symmetry::Name& newSymmetry,
  const ChiralStatePreservation& preservationOption
) {
  auto canonicalizedSubstituents = glue::canonicalize(newRanking.sortedSubstituents);
  auto canonicalizedCharacters = glue::makeCanonicalCharacters(canonicalizedSubstituents);
  auto newLinks = glue::makeLinksSelfReferential(canonicalizedSubstituents, newRanking.links);
  auto newStereopermutations = stereopermutation::uniqueStereopermutationsWithWeights(
    StereopermutationType {
      newSymmetry,
      canonicalizedCharacters,
      newLinks
    },
    newSymmetry
  );

  _removeImpossibleStereopermutations(newStereopermutations, molecule);

  boost::optional<unsigned> newStereopermutation;

  /* Does the stereocenter carry chiral information?
   * If so, try to get a symmetry mapping to the new symmetry position
   * If there are mappings, try to select one according to preservationOption policy
   *
   * If any of those steps returns boost::none, the whole expression is
   * boost::none.
   */
  auto suitableMappingOptional = _assignmentOption
    | temple::callIfSome(
        Symmetry::getMapping,
        _symmetry,
        newSymmetry,
        // Last parameter is the deleted symmetry position, get this from cache
        _symmetryPositionMapCache.at(which)
      )
    | temple::callIfSome(
        glue::getIndexMapping,
        temple::ANS,
        preservationOption
      );

  if(suitableMappingOptional) {
    const auto& symmetryMapping = suitableMappingOptional.value();

    // Transform current assignment from characters to atom indices
    const auto& currentStereopermutation = _uniqueStereopermutationsCache.assignments.at(
      _assignmentOption.value()
    );

    std::vector<AtomIndexType> atomsAtCurrentSymmetryPositions = glue::mapToSymmetryPositions(
      currentStereopermutation,
      _ranking.sortedSubstituents
    );

    // Transfer indices from current symmetry to new symmetry
    std::vector<AtomIndexType> atomsAtNewSymmetryPositions (Symmetry::size(newSymmetry));
    for(unsigned i = 0; i < Symmetry::size(newSymmetry); ++i) {
      atomsAtNewSymmetryPositions.at(
        symmetryMapping.at(i)
      ) = atomsAtCurrentSymmetryPositions.at(i);
    }

    // Get character representation in new symmetry
    std::vector<char> charactersInNewSymmetry = glue::makeStereopermutationCharacters(
      canonicalizedSubstituents,
      canonicalizedCharacters,
      atomsAtNewSymmetryPositions
    );

    // Construct an assignment
    auto trialStereopermutation = stereopermutation::Stereopermutation(
      newSymmetry,
      charactersInNewSymmetry,
      newLinks
    );

    // Generate the rotational equivalents
    auto allTrialRotations = trialStereopermutation.generateAllRotations(newSymmetry);

    // Search for a match from the vector of uniques
    for(unsigned i = 0; i < newStereopermutations.assignments.size(); ++i) {
      if(allTrialRotations.count(newStereopermutations.assignments.at(i)) > 0) {
        newStereopermutation = i;
        break;
      }
    }
  }

  // Overwrite class state
  _ranking = newRanking;
  _ranking.sortedSubstituents = canonicalizedSubstituents;
  _symmetry = newSymmetry;
  _uniqueStereopermutationsCache = newStereopermutations;
  assign(newStereopermutation);
}

Symmetry::Name CNStereocenter::getSymmetry() const {
  return _symmetry;
}

void CNStereocenter::fit(
  const Molecule& molecule,
  const Delib::PositionCollection& positions,
  std::vector<Symmetry::Name> excludeSymmetries
) {
  // Extract a list of adjacent indices from stored state
  std::vector<AtomIndexType> adjacentAtoms;

  for(const auto& equalPrioritySet : _ranking.sortedSubstituents) {
    for(const auto& substituentIndex : equalPrioritySet) {
      adjacentAtoms.push_back(substituentIndex);
    }
  }

  const Symmetry::Name priorSymmetry {this->_symmetry};
  const boost::optional<unsigned> priorStereopermutation {this->_assignmentOption};

  const Symmetry::Name initialSymmetry {Symmetry::Name::Linear};
  const unsigned initialStereopermutation = 0;
  const unsigned initialPenalty = 100;

  Symmetry::Name bestSymmetry = initialSymmetry;
  unsigned bestStereopermutation = initialStereopermutation;
  double bestPenalty = initialPenalty;
  unsigned bestStereopermutationMultiplicity = 1;

  auto excludesContains = temple::makeContainsPredicate(excludeSymmetries);

  // Cycle through all symmetries
  for(const auto& symmetryName : Symmetry::allNames) {
    // Skip any Symmetries of different size
    if(
      Symmetry::size(symmetryName) != Symmetry::size(_symmetry)
      || excludesContains(symmetryName)
    ) {
      continue;
    }

    // Change the symmetry of the CNStereocenter
    setSymmetry(molecule, symmetryName);

    for(
      unsigned assignment = 0;
      assignment < numStereopermutations();
      ++assignment
    ) {
      // Assign the stereocenter
      assign(assignment);

      const auto prototypes = chiralityConstraints();

      const double angleDeviations = temple::sum(
        temple::mapAllPairs(
          adjacentAtoms,
          [&](const AtomIndexType& i, const AtomIndexType& k) -> double {
            return std::fabs(
              DelibHelpers::getAngle( // The angle from the positions
                positions,
                i,
                _centerAtom,
                k
              ) - angle(
                i,
                _centerAtom,
                k
              )
            );
          }
        )
      );

      // We can stop immediately if this is worse
      if(angleDeviations > bestPenalty) {
        continue;
      }

      const double oneThreeDistanceDeviations = temple::sum(
        temple::mapAllPairs(
          adjacentAtoms,
          [&](const AtomIndexType& i, const AtomIndexType& k) -> double {
            return std::fabs(
              DelibHelpers::getDistance( // i-k 1-3 distance from positions
                positions,
                i,
                k
              ) - CommonTrig::lawOfCosines( // idealized 1-3 distance from
                DelibHelpers::getDistance( // i-j 1-2 distance from positions
                  positions,
                  i,
                  _centerAtom
                ),
                DelibHelpers::getDistance( // j-k 1-2 distance from positions
                  positions,
                  _centerAtom,
                  k
                ),
                angle( // idealized Stereocenter angle
                  i,
                  _centerAtom,
                  k
                )
              )
            );
          }
        )
      );

      // Another early continue
      if(angleDeviations + oneThreeDistanceDeviations > bestPenalty) {
        continue;
      }

      const double chiralityDeviations = (prototypes.empty()
        ? 0
        : temple::sum(
          temple::map(
            prototypes,
            [&positions](const auto& constraintPrototype) -> double {
              using TargetEnumType = Stereocenters::ChiralityConstraintTarget;

              double volume = DelibHelpers::getSignedVolume(
                positions,
                constraintPrototype.indices
              );

              // If the target is flat, then the "error" is continuous:
              if(constraintPrototype.target == TargetEnumType::Flat) {
                return std::fabs(volume);
              }

              /* Otherwise, no bounds, error is some arbitrary penalty if the sign is
               * wrong
               */
              if(
                (
                  constraintPrototype.target == TargetEnumType::Positive
                  && volume < 0
                ) || (
                  constraintPrototype.target == TargetEnumType::Negative
                  && volume > 0
                )
              ) {
                return 1; // Arbitrary penalty
              }

              return 0;
            }
          )
        )
      );

      double fitPenalty = angleDeviations
        + oneThreeDistanceDeviations
        + chiralityDeviations;


#ifndef NDEBUG
      Log::log(Log::Particulars::CNStereocenterFit)
        << Symmetry::nameIndex(symmetryName)
        << ", " << assignment
        << ", " << std::setprecision(4) << std::fixed
        << angleDeviations << ", "
        << oneThreeDistanceDeviations << ", "
        << chiralityDeviations
        << std::endl;
#endif

      if(fitPenalty < bestPenalty) {
        bestSymmetry = symmetryName;
        bestStereopermutation = assignment;
        bestPenalty = fitPenalty;
        bestStereopermutationMultiplicity = 1;
      } else if(fitPenalty == bestPenalty) {
        // Assume that IF we have multiplicity, it's from the same symmetry
        assert(bestSymmetry == symmetryName);
        bestStereopermutationMultiplicity += 1;
      }
    }
  }

  /* In case NO assignments could be tested, return to the prior state.
   * This guards against situations in which predicates in uniqueStereopermutations
   * could lead no assignments to be returned, such as in e.g. square-planar
   * AAAB with {0, 3}, {1, 3}, {2, 3} with removal of trans-spanning groups.
   * In that situation, all possible assignments are trans-spanning and
   * uniqueStereopermutations is an empty vector.
   *
   * At the moment, this predicate is disabled, so no such issues should arise.
   * Just being safe.
   */
  if(
    bestSymmetry == initialSymmetry
    && bestStereopermutation == initialStereopermutation
    && bestPenalty == initialPenalty
  ) {
    // Return to prior
    setSymmetry(molecule, priorSymmetry);
    assign(priorStereopermutation);

  } else {
    // Set to best fit
    setSymmetry(molecule, bestSymmetry);

    /* How to handle multiplicity?
     * Current policy: If there is multiplicity, warn and do not assign
     */
    if(bestStereopermutationMultiplicity > 1) {
      assign(boost::none);
    } else {
      assign(bestStereopermutation);
    }
  }
}

/* Information */
double CNStereocenter::angle(
  const AtomIndexType& i,
  const AtomIndexType& j __attribute__((unused)),
  const AtomIndexType& k
) const {
  assert(j == _centerAtom);

  /* j is practically unused here because in the Symmetry angleFunctions, the
   * middle atom is implicit, it has no symmetryPosition number. j is however
   * needed in the interface because of EZStereocenter, where it is important
   * to specify which atom is the intermediate (there are two possibilities)
   */
  return Symmetry::angleFunction(_symmetry)(
    _symmetryPositionMapCache.at(i),
    _symmetryPositionMapCache.at(k)
  );
}

#if false
DistanceGeometry::ValueBounds CNStereocenter::angles(
  const AtomIndexType& i,
  const AtomIndexType& j __attribute__((unused)),
  const AtomIndexType& k
) const {
  assert(j == _centerAtom);

  // Require:
  // Mapping from atom index to ligand index
  std::map<AtomIndexType, unsigned> _ligandMap;

  /* j is practically unused here because in the Symmetry angleFunctions, the
   * middle atom is implicit, it has no symmetryPosition number. j is however
   * needed in the interface because of EZStereocenter, where it is important
   * to specify which atom is the intermediate (there are two possibilities)
   */
  assert(_ligandMap.at(i) != _ligandMap.at(k));


}
#endif


boost::optional<unsigned> CNStereocenter::assigned() const {
  return _assignmentOption;
}

std::vector<
  ChiralityConstraintPrototype
> CNStereocenter::chiralityConstraints() const {
  std::vector<ChiralityConstraintPrototype> prototypes;

  // Only collect constraints if there is more than one assignment for this
  if(numStereopermutations() > 1) {

    /* Invert _neighborSymmetryPositionMap, we need a mapping of
     *  (position in symmetry) -> atom index
     */
    auto symmetryPositionToAtomIndexMap = temple::invertMap(
      _symmetryPositionMapCache
    );

    // Get list of tetrahedra from symmetry
    auto tetrahedraList = Symmetry::tetrahedra(_symmetry);

    for(const auto& tetrahedron : tetrahedraList) {
      /* Replace boost::none with centerAtom, indices (represent positions within
       * the symmetry) with the atom index at that position from the inverted map
       */
      auto replaced = temple::map(
        tetrahedron,
        [&](const auto& indexOptional) -> AtomIndexType {
          if(indexOptional) {
            return symmetryPositionToAtomIndexMap.at(indexOptional.value());
          }

          return _centerAtom;
        }
      );

      // Make a prototype from it
      prototypes.emplace_back(
        replaced,
        ChiralityConstraintTarget::Positive
      );
    }

  }

  return prototypes;
}

std::vector<DihedralLimits> CNStereocenter::dihedralLimits() const {
  return {};
}

std::string CNStereocenter::info() const {
  // TODO revisit as soon as linking information is introduced
  std::string returnString = "CN "s
    + std::to_string(_centerAtom) + " ("s + Symmetry::name(_symmetry) +", "s;

  auto characters = glue::makeCanonicalCharacters(_ranking.sortedSubstituents);
  std::copy(
    characters.begin(),
    characters.end(),
    std::back_inserter(returnString)
  );

  for(const auto& link : glue::makeLinksSelfReferential(_ranking.sortedSubstituents, _ranking.links)) {
    returnString += ", "s + characters.at(link.first) + "-"s + characters.at(link.second);
  }

  returnString += "): "s;

  if(_assignmentOption) {
    returnString += std::to_string(_assignmentOption.value());
  } else {
    returnString += "u";
  }

  returnString += "/"s + std::to_string(numStereopermutations());

  return returnString;
}

std::string CNStereocenter::rankInfo() const {
  return (
    "CN-"s + std::to_string(static_cast<unsigned>(_symmetry))
    + "-"s + std::to_string(numStereopermutations())
    + "-"s + (
      assigned()
      ? std::to_string(assigned().value())
      : "u"s
    )
  );
}

std::vector<AtomIndexType> CNStereocenter::involvedAtoms() const {
  return {_centerAtom};
}

unsigned CNStereocenter::numStereopermutations() const {
  return _uniqueStereopermutationsCache.assignments.size();
}

void CNStereocenter::setSymmetry(
  const Molecule& molecule,
  const Symmetry::Name& symmetryName
) {
  // Set new symmetry
  _symmetry = symmetryName;

  // recalculate the number of unique Stereopermutations
  _uniqueStereopermutationsCache = stereopermutation::uniqueStereopermutationsWithWeights(
    StereopermutationType(
      symmetryName,
      glue::makeCanonicalCharacters(_ranking.sortedSubstituents),
      glue::makeLinksSelfReferential(_ranking.sortedSubstituents, _ranking.links)
    ),
    symmetryName,
    false // do NOT remove trans-spanning ligand groups
  );

  _removeImpossibleStereopermutations(_uniqueStereopermutationsCache, molecule);

  _symmetryPositionMapCache.clear();

  // The Stereocenter is now unassigned
  _assignmentOption = boost::none;
}

Type CNStereocenter::type() const {
  return Type::CNStereocenter;
}

bool CNStereocenter::operator == (const CNStereocenter& other) const {
  return (
    _symmetry == other._symmetry
    && _centerAtom == other._centerAtom
    && _uniqueStereopermutationsCache.assignments.size() == other._uniqueStereopermutationsCache.assignments.size()
    && _assignmentOption == other._assignmentOption
  );
}

bool CNStereocenter::operator < (const CNStereocenter& other) const {
  /* Sequentially compare individual components, comparing assignments last
   * if everything else matches
   */
  return temple::consecutiveCompareSmaller(
    _centerAtom,
    other._centerAtom,
    _uniqueStereopermutationsCache.assignments.size(),
    other._uniqueStereopermutationsCache.assignments.size(),
    _symmetry,
    other._symmetry,
    _assignmentOption,
    other._assignmentOption
  );
}

} // namespace Stereocenters

} // namespace molassembler
