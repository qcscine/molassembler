#include "template_magic/Boost.h"
#include "template_magic/Containers.h"
#include "template_magic/Numeric.h"
#include "template_magic/Random.h"
#include "geometry_assignment/GenerateUniques.h"
#include "symmetry_information/Properties.h"

#include "BuildTypeSwitch.h"
#include "CNStereocenter.h"
#include "CommonTrig.h"
#include "DelibHelpers.h"
#include "Log.h"
#include "StdlibTypeAlgorithms.h"

#include <iomanip>

/* TODO
 * - An option to enable non-minimal-work mappings if unique
 * - Lots of parallels between addSubstituent, propagateGraphChange and
 *   removeSubstituent. Refactor? Very difficult, at least a ten-parameter
 *   function with optionals to smooth over differences...
 */

namespace MoleculeManip {

namespace Stereocenters {

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
  const RankingInformation::RankedType canonicalizedSubstituents
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

UniqueAssignments::Assignment::LinksSetType makeLinksSelfReferential(
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

  UniqueAssignments::Assignment::LinksSetType links;

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
    links.emplace(
      findIndexInSorted(linkPair.first),
      findIndexInSorted(linkPair.second)
    );
  }

  return links;
}

std::map<AtomIndexType, unsigned> makeSymmetryPositionMap(
  const UniqueAssignments::Assignment& assignment,
  const RankingInformation::RankedType& canonicalizedSubstituents
) {
  std::map<AtomIndexType, unsigned> positionMap;

  const auto endIterators = TemplateMagic::map(
    canonicalizedSubstituents,
    [](const auto& equalPrioritySet) -> auto {
      return equalPrioritySet.end();
    }
  );

  auto iterators = TemplateMagic::map(
    canonicalizedSubstituents,
    [](const auto& equalPrioritySet) -> auto {
      return equalPrioritySet.begin();
    }
  );

  unsigned symmetryPosition = 0;
  for(const auto& priorityChar : assignment.characters) {
    auto& rankingIterator = iterators.at(priorityChar - 'A');
    const auto& endIterator = endIterators.at(priorityChar - 'A');
    assert(rankingIterator != endIterator);

    positionMap.emplace(
      *rankingIterator,
      symmetryPosition
    );

    ++symmetryPosition;
    ++rankingIterator;
  }

  return positionMap;
}


std::vector<AtomIndexType> mapToSymmetryPositions(
  const UniqueAssignments::Assignment& assignment,
  const RankingInformation::RankedType& canonicalizedSubstituents
) {
  std::vector<AtomIndexType> atomsAtSymmetryPositions;

  const auto endIterators = TemplateMagic::map(
    canonicalizedSubstituents,
    [](const auto& equalPrioritySet) {
      return equalPrioritySet.end();
    }
  );

  auto iterators = TemplateMagic::map(
    canonicalizedSubstituents,
    [](const auto& equalPrioritySet) -> auto {
      return equalPrioritySet.begin();
    }
  );

  for(const auto& priorityChar : assignment.characters) {
    auto& rankingIterator = iterators.at(priorityChar - 'A');
    const auto& endIterator = endIterators.at(priorityChar - 'A');
    assert(rankingIterator != endIterator);

    atomsAtSymmetryPositions.push_back(
      *rankingIterator
    );

    ++rankingIterator;
  }

  return atomsAtSymmetryPositions;
}

std::vector<char> makeAssignmentCharacters(
  const RankingInformation::RankedType& canonicalizedSubstituents,
  const std::vector<char>& canonicalizedAssignmentCharacters,
  const std::vector<AtomIndexType>& atomsAtSymmetryPositions
) {
  // Replace the atom indices by their new ranking characters
  std::vector<unsigned> flattenedIndices;
  for(const auto& equalPrioritySet : canonicalizedSubstituents) {
    for(const auto& index : equalPrioritySet) {
      flattenedIndices.push_back(index);
    }
  }

  std::vector<char> newAssignmentCharacters;

  for(const auto& index : atomsAtSymmetryPositions) {
    const auto findIter = std::find(
      flattenedIndices.begin(),
      flattenedIndices.end(),
      index
    );

    newAssignmentCharacters.push_back(
      canonicalizedAssignmentCharacters.at(
        findIter - flattenedIndices.begin()
      )
    );
  }

  return newAssignmentCharacters;
}

boost::optional<const std::vector<unsigned>&> getIndexMapping(
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
      TemplateMagic::random.getSingle<unsigned>(
        0,
        mappingsGroup.indexMappings.size() - 1
      )
    );
  }

  return boost::none;
}

} // namespace glue

/* Constructors */
CNStereocenter::CNStereocenter(
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
  _uniqueAssignmentsCache = UniqueAssignments::uniqueAssignments(
    AssignmentType(
      symmetry,
      glue::makeCanonicalCharacters(_ranking.sortedSubstituents),
      glue::makeLinksSelfReferential(_ranking.sortedSubstituents, _ranking.linkedPairs)
    ),
    false // Do NOT remove trans-spanning linked groups
  );
}

/* Modification */
void CNStereocenter::addSubstituent(
  const AtomIndexType& newSubstituentIndex,
  const RankingInformation& newRanking,
  const Symmetry::Name& newSymmetry,
  const ChiralStatePreservation& preservationOption
) {
  auto canonicalizedSubstituents = glue::canonicalize(newRanking.sortedSubstituents);
  auto canonicalizedCharacters = glue::makeCanonicalCharacters(canonicalizedSubstituents);
  auto newLinks = glue::makeLinksSelfReferential(canonicalizedSubstituents, newRanking.linkedPairs);
  auto newAssignments = UniqueAssignments::uniqueAssignments(
    AssignmentType {
      newSymmetry,
      canonicalizedCharacters,
      newLinks
    }
  );
  boost::optional<unsigned> newAssignment = boost::none;

  // Does the current stereocenter carry chiral information?
  if(_assignmentOption) {
    auto symmetryMappingOptional = Symmetry::getMapping(
      _symmetry,
      newSymmetry
    );
    // Are there mappings to the target symmetry?
    if(symmetryMappingOptional) {
      auto suitableMappingOptional = glue::getIndexMapping(
        symmetryMappingOptional.value(),
        preservationOption
      );
      // Is there a suitable mapping?
      if(suitableMappingOptional) {
        /* So now we must transfer the current assignment into the new symmetry
         * and search for it in the set of uniques.
         */
        const auto& symmetryMapping = suitableMappingOptional.value();

        // Transform current assignment from characters to indices
        const auto& currentAssignment = _uniqueAssignmentsCache.at(
          _assignmentOption.value()
        );

        std::vector<AtomIndexType> atomsAtSmallerSymmetryPositions = glue::mapToSymmetryPositions(
          currentAssignment,
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
        std::vector<char> charactersInLargerSymmetry = glue::makeAssignmentCharacters(
          canonicalizedSubstituents,
          canonicalizedCharacters,
          atomsAtLargerSymmetryPositions
        );

        // Construct an assignment from it
        auto trialAssignment = UniqueAssignments::Assignment(
          newSymmetry,
          charactersInLargerSymmetry,
          newLinks
        );

        // Generate the rotational equivalents
        auto allTrialRotations = trialAssignment.generateAllRotations();

        // Search for a match from the vector of uniques
        for(unsigned i = 0; i < newAssignments.size(); ++i) {
          if(allTrialRotations.count(newAssignments.at(i)) > 0) {
            newAssignment = i;
            break;
          }
        }
      }
    };
  } 

  /* Since either there is no steric information present or no unique
   * minimal-effort mapping exists to the new symmetry, it is impossible to
   * unambiguously choose a new assignment. In those cases, newAssignment
   * remains boost::none, and this stereocenter loses any chiral information it
   * may have had.
   */

  // Overwrite class state
  _ranking = newRanking;
  _ranking.sortedSubstituents = canonicalizedSubstituents;
  _symmetry = newSymmetry;
  _uniqueAssignmentsCache = newAssignments;
  assign(newAssignment);
}

void CNStereocenter::assign(const boost::optional<unsigned>& assignment) {
  if(assignment) {
    assert(assignment.value() < _uniqueAssignmentsCache.size());
  }

  // Store current assignment
  _assignmentOption = assignment;

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndexType -> unsigned).
   */
  if(assignment) {
    _symmetryPositionMapCache = glue::makeSymmetryPositionMap(
      _uniqueAssignmentsCache.at(assignment.value()),
      _ranking.sortedSubstituents
    );
  } else { // Wipe the map
    _symmetryPositionMapCache.clear();
  }
}

void CNStereocenter::propagateGraphChange(const RankingInformation& newRanking) {
  auto newCanonicalizedSubstituents = glue::canonicalize(newRanking.sortedSubstituents);

  // Has anything about the new ranking changed relative to the old?
  if(
    newCanonicalizedSubstituents != _ranking.sortedSubstituents
    || newRanking.linkedPairs != _ranking.linkedPairs
  ) {
    auto newCharacters = glue::makeCanonicalCharacters(newCanonicalizedSubstituents);
    auto newLinks = glue::makeLinksSelfReferential(
      newCanonicalizedSubstituents,
      newRanking.linkedPairs
    );
    auto newAssignments = UniqueAssignments::uniqueAssignments(
      AssignmentType(_symmetry, newCharacters, newLinks),
      false // do NOT remove trans-spanning ligand groups
    );
    boost::optional<unsigned> newAssignment = boost::none;

    /* Before we overwrite class state, we need to figure out which assignment 
     * in the new set of assignments corresponds to the one we have now.
     * This is only necessary in the case that the stereocenter is currently
     * assigned and only possible if the new number of assignments is smaller or
     * equal to the amount we have currently.
     */
    if(
      _assignmentOption
      && newAssignments.size() <= _uniqueAssignmentsCache.size()
    ) {
      const auto& currentAssignment = _uniqueAssignmentsCache.at(
        _assignmentOption.value()
      );

      // Replace the characters by their corresponding indices from the old ranking
      std::vector<AtomIndexType> atomsAtSymmetryPositions = glue::mapToSymmetryPositions(
        currentAssignment,
        _ranking.sortedSubstituents
      );

      // Replace the atom indices by their new ranking characters
      std::vector<char> newAssignmentCharacters = glue::makeAssignmentCharacters(
        newCanonicalizedSubstituents,
        newCharacters,
        atomsAtSymmetryPositions
      );

      // Create a new assignment with those characters
      auto trialAssignment = UniqueAssignments::Assignment(
        _symmetry,
        newAssignmentCharacters,
        newLinks
      );

      // Generate all rotations of this trial assignment
      auto allTrialRotations = trialAssignment.generateAllRotations();

      // Find out which of the new assignments has a rotational equivalent
      for(unsigned i = 0; i < newAssignments.size(); ++i) {
        if(allTrialRotations.count(newAssignments.at(i)) > 0) {
          newAssignment = i;
          break;
        }
      }
    }

    // Overwrite the class state
    _ranking = newRanking;
    _ranking.sortedSubstituents = newCanonicalizedSubstituents;
    _uniqueAssignmentsCache = newAssignments;
    assign(newAssignment);
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

  RankingInformation::LinksType newLinks;
  for(auto& linkedPair : _ranking.linkedPairs) {
    newLinks.emplace(
      updateIndex(linkedPair.first),
      updateIndex(linkedPair.second)
    );
  }
  _ranking.linkedPairs = std::move(newLinks);

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
  const AtomIndexType& which,
  const RankingInformation& newRanking,
  const Symmetry::Name& newSymmetry,
  const ChiralStatePreservation& preservationOption
) {
  auto canonicalizedSubstituents = glue::canonicalize(newRanking.sortedSubstituents);
  auto canonicalizedCharacters = glue::makeCanonicalCharacters(canonicalizedSubstituents);
  auto newLinks = glue::makeLinksSelfReferential(canonicalizedSubstituents, newRanking.linkedPairs);
  auto newAssignments = UniqueAssignments::uniqueAssignments(
    AssignmentType {
      newSymmetry,
      canonicalizedCharacters,
      newLinks
    }
  );
  boost::optional<unsigned> newAssignment;


  // Does the stereocenter currently carry chiral information?
  if(_assignmentOption) {
    // Which symmetry position is the deleted index currently at?
    unsigned deletedSymmetryPosition = _symmetryPositionMapCache.at(which);

    // Fetch the mapping group optional
    auto symmetryMappingOptional = Symmetry::getMapping(
      _symmetry,
      newSymmetry,
      deletedSymmetryPosition
    );

    if(symmetryMappingOptional) {
      auto suitableMappingOptional = glue::getIndexMapping(
        symmetryMappingOptional.value(),
        preservationOption
      );

      if(suitableMappingOptional) {
        const auto& symmetryMapping = suitableMappingOptional.value();

        // Transform current assignment from characters to atom indices
        const auto& currentAssignment = _uniqueAssignmentsCache.at(
          _assignmentOption.value()
        );

        std::vector<AtomIndexType> atomsAtCurrentSymmetryPositions = glue::mapToSymmetryPositions(
          currentAssignment,
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
        std::vector<char> charactersInNewSymmetry = glue::makeAssignmentCharacters(
          canonicalizedSubstituents,
          canonicalizedCharacters,
          atomsAtNewSymmetryPositions
        );

        // Construct an assignment
        auto trialAssignment = UniqueAssignments::Assignment(
          newSymmetry,
          charactersInNewSymmetry,
          newLinks
        );

        // Generate the rotational equivalents
        auto allTrialRotations = trialAssignment.generateAllRotations();

        // Search for a match from the vector of uniques
        for(unsigned i = 0; i < newAssignments.size(); ++i) {
          if(allTrialRotations.count(newAssignments.at(i)) > 0) {
            newAssignment = i;
            break;
          }
        }
      }
    }
  }

  // Overwrite class state
  _ranking = newRanking;
  _ranking.sortedSubstituents = canonicalizedSubstituents;
  _symmetry = newSymmetry;
  _uniqueAssignmentsCache = newAssignments;
  assign(newAssignment);
}

const AtomIndexType& CNStereocenter::getCentralAtomIndex() const {
  return _centerAtom;
}

Symmetry::Name CNStereocenter::getSymmetry() const {
  return _symmetry;
}

void CNStereocenter::fit(
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
  const boost::optional<unsigned> priorAssignment {this->_assignmentOption};

  const Symmetry::Name initialSymmetry {Symmetry::Name::Linear};
  const unsigned initialAssignment = 0;
  const unsigned initialPenalty = 100;

  Symmetry::Name bestSymmetry = initialSymmetry;
  unsigned bestAssignment = initialAssignment;
  double bestPenalty = initialPenalty;
  unsigned bestAssignmentMultiplicity = 1;

  auto excludesContains = TemplateMagic::makeContainsPredicate(excludeSymmetries);

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
    setSymmetry(symmetryName);

    for(
      unsigned assignment = 0;
      assignment < numAssignments();
      assignment++
    ) {
      // Assign the stereocenter
      assign(assignment);

      const auto prototypes = chiralityConstraints();

      const double angleDeviations = TemplateMagic::sum(
        TemplateMagic::mapAllPairs(
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

      const double oneThreeDistanceDeviations = TemplateMagic::sum(
        TemplateMagic::mapAllPairs(
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

      const double chiralityDeviations = (prototypes.empty()
        ? 0
        : TemplateMagic::sum(
          TemplateMagic::map(
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


      /* If ever you want to attempt to additionally penalize symmetries that
       * the internal determineSymmetry functionality disfavors, then that
       * information must come from a different prepending function call, it is
       * not a valid part of the Stereocenter interface... EZStereocenter
       * cannot expect to be penalized in the same way
       */
      /* if(
        expectedSymmetry 
        && expectedSymmetry.value() != symmetryName
      ) {
        currentFit.symmetryPenalty = 0.5;
      } */

      if(fitPenalty < bestPenalty) {
        bestSymmetry = symmetryName;
        bestAssignment = assignment;
        bestPenalty = fitPenalty;
        bestAssignmentMultiplicity = 1;
      } else if(fitPenalty == bestPenalty) {
        // Assume that IF we have multiplicity, it's from the same symmetry
        assert(bestSymmetry == symmetryName);
        bestAssignmentMultiplicity += 1;
      }
    }
  }
  
  /* In case NO assignments could be tested, return to the prior state.
   * This guards against situations in which predicates in uniqueAssignments
   * could lead no assignments to be returned, such as in e.g. square-planar
   * AAAB with {0, 3}, {1, 3}, {2, 3} with removal of trans-spanning groups.
   * In that situation, all possible assignments are trans-spanning and 
   * uniqueAssignments is an empty vector.
   *
   * At the moment, this predicate is disabled, so no such issues should arise.
   * Just being safe.
   */
  if( 
    bestSymmetry == initialSymmetry
    && bestAssignment == initialAssignment
    && bestPenalty == initialPenalty
  ) {
    // Return to prior
    setSymmetry(priorSymmetry);
    assign(priorAssignment);

  } else {
    // Set to best fit
    setSymmetry(bestSymmetry);

    /* How to handle multiplicity? 
     * Current policy: If there is multiplicity, warn and do not assign
     */
    if(bestAssignmentMultiplicity > 1) {
      assign(boost::none);
    } else {
      assign(bestAssignment);
    }
  }
}

/* Information */
double CNStereocenter::angle(
  const AtomIndexType& i,
  const AtomIndexType& j, 
  const AtomIndexType& k
) const {
  assert(j == _centerAtom);

  /* j is pracitcally unused here because in the Symmetry angleFunctions, the
   * middle atom is implicit, it has no symmetryPosition number. j is however
   * needed in the interface because of EZStereocenter, where it is important
   * to specify which atom is the intermediate (there are two possibilities)
   */
  return Symmetry::angleFunction(_symmetry)(
    _symmetryPositionMapCache.at(i),
    _symmetryPositionMapCache.at(k)
  );
}

boost::optional<unsigned> CNStereocenter::assigned() const {
  return _assignmentOption;
}

std::vector<
  ChiralityConstraintPrototype
> CNStereocenter::chiralityConstraints() const {
  std::vector<ChiralityConstraintPrototype> prototypes;

  // Only collect constraints if there is more than one assignment for this
  if(numAssignments() > 1) {
  
    /* Invert _neighborSymmetryPositionMap, we need a mapping of
     *  (position in symmetry) -> atom index
     */
    auto symmetryPositionToAtomIndexMap = TemplateMagic::invertMap(
      _symmetryPositionMapCache
    );

    // Get list of tetrahedra from symmetry
    auto tetrahedraList = Symmetry::tetrahedra(_symmetry);

    for(const auto& tetrahedron : tetrahedraList) {
      /* Replace boost::none with centerAtom, indices (represent positions within 
       * the symmetry) with the atom index at that position from the inverted map
       */
      auto replaced = TemplateMagic::map(
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

  returnString += "): "s;

  if(_assignmentOption) {
    returnString += std::to_string(_assignmentOption.value());
  } else {
    returnString += "u";
  }

  returnString += "/"s + std::to_string(numAssignments());

  return returnString;
}

std::string CNStereocenter::rankInfo() const {
  // TODO revisit as soon as pseudo-asymmetry is added

  return (
    "CN-"s + std::to_string(static_cast<unsigned>(_symmetry)) 
    + "-"s + std::to_string(numAssignments())
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

unsigned CNStereocenter::numAssignments() const {
  return _uniqueAssignmentsCache.size();
}

void CNStereocenter::setSymmetry(const Symmetry::Name& symmetryName) {
  // Set new symmetry
  _symmetry = symmetryName;

  // recalculate the number of unique Assignments
  _uniqueAssignmentsCache = UniqueAssignments::uniqueAssignments(
    AssignmentType(
      symmetryName,
      glue::makeCanonicalCharacters(_ranking.sortedSubstituents),
      glue::makeLinksSelfReferential(_ranking.sortedSubstituents, _ranking.linkedPairs)
    ),
    false // do NOT remove trans-spanning ligand groups
  );

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
    && _uniqueAssignmentsCache.size() == other._uniqueAssignmentsCache.size()
    && _assignmentOption == other._assignmentOption
  );
}

bool CNStereocenter::operator < (const CNStereocenter& other) const {
  using TemplateMagic::componentSmaller;
  
  /* Sequentially compare individual components, comparing assignments last
   * if everything else matches
   */
  return componentSmaller(
    _centerAtom,
    other._centerAtom
  ).value_or(
    componentSmaller(
      _uniqueAssignmentsCache.size(),
      other._uniqueAssignmentsCache.size()
    ).value_or(
      componentSmaller(
        _symmetry,
        other._symmetry
      ).value_or(
        // NOTE: boost::none is smaller than 0 in this mixed ordering
        _assignmentOption < other._assignmentOption
      )
    )
  );
}

} // namespace Stereocenters

} // namespace MoleculeManip
