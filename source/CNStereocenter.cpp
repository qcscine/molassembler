#include "template_magic/Boost.h"
#include "template_magic/Containers.h"
#include "template_magic/Numeric.h"
#include "geometry_assignment/GenerateUniques.h"

#include "CNStereocenter.h"
#include "CommonTrig.h"
#include "DelibHelpers.h"
#include "Log.h"
#include "StdlibTypeAlgorithms.h"

#include <iomanip>

namespace MoleculeManip {

namespace Stereocenters {

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
  // Generate the set of unique assignments possible here
  _uniqueAssignmentsCache = UniqueAssignments::uniqueAssignments(
    AssignmentType(
      symmetry,
      _makeAssignmentCharacters(_ranking),
      _makeAssignmentLinks(_ranking)
    ),
    false // Do NOT remove trans-spanning linked groups
  );
}

/* Private members */
std::vector<char> CNStereocenter::_makeAssignmentCharacters(
  const RankingInformation& ranking
) {
  std::vector<char> characters;

  char currentChar = 'A';
  for(const auto& equalPrioritySet : ranking.sortedSubstituents) {
    for(unsigned i = 0; i < equalPrioritySet.size(); ++i) {
      characters.push_back(currentChar);
    }

    ++currentChar;
  }

  return characters;
}

std::map<AtomIndexType, unsigned> CNStereocenter::_makeSymmetryPositionMap(
  const UniqueAssignments::Assignment& assignment,
  const RankingInformation& ranking
) {
  std::map<AtomIndexType, unsigned> positionMap;

  const auto endIterators = TemplateMagic::map(
    ranking.sortedSubstituents,
    [](const auto& equalPrioritySet) -> auto {
      return equalPrioritySet.end();
    }
  );

  auto iterators = TemplateMagic::map(
    ranking.sortedSubstituents,
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

UniqueAssignments::Assignment::LinksSetType CNStereocenter::_makeAssignmentLinks(
  const RankingInformation& ranking
) {
  // Flatten the sorted list of indices
  std::vector<unsigned> sortedIndices;
  for(const auto& equalPrioritySet : ranking.sortedSubstituents) {
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

  for(const auto& linkPair : ranking.linkedPairs) {
    links.emplace(
      findIndexInSorted(linkPair.first),
      findIndexInSorted(linkPair.second)
    );
  }

  return links;
}

/* Modification */
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
    _symmetryPositionMapCache = _makeSymmetryPositionMap(
      _uniqueAssignmentsCache.at(assignment.value()),
      _ranking
    );
  } else { // Wipe the map
    _symmetryPositionMapCache.clear();
  }
}

void CNStereocenter::adaptToRankingChange(const RankingInformation& newRanking) {
  /* TODO previous symmetry is NOT taken into account -> no mapping to new
   * symmetry possible
   *
   * This is a mess, and there is no clear logic to it
   *
   * This is for when something is EDITED in an existing molecule, this has no
   * priority
   */

}

const AtomIndexType& CNStereocenter::getCentralAtomIndex() const {
  return _centerAtom;
}

Symmetry::Name CNStereocenter::getSymmetry() const {
  return _symmetry;
}

void CNStereocenter::changeSymmetry(const Symmetry::Name& symmetryName) {
  // Set new symmetry
  _symmetry = symmetryName;

  // recalculate the number of unique Assignments
  _uniqueAssignmentsCache = UniqueAssignments::uniqueAssignments(
    AssignmentType(
      symmetryName,
      _makeAssignmentCharacters(_ranking),
      _makeAssignmentLinks(_ranking)
    ),
    false // do NOT remove trans-spanning ligand groups
  );

  _symmetryPositionMapCache.clear();

  // The Stereocenter is now unassigned
  _assignmentOption = boost::none;
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
    changeSymmetry(symmetryName);

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
    changeSymmetry(priorSymmetry);
    assign(priorAssignment);

  } else {
    // Set to best fit
    changeSymmetry(bestSymmetry);

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

  auto characters = _makeAssignmentCharacters(_ranking);
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

std::set<AtomIndexType> CNStereocenter::involvedAtoms() const {
  return {_centerAtom};
}

unsigned CNStereocenter::numAssignments() const {
  return _uniqueAssignmentsCache.size();
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
