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
) : _neighborCharMap {
      _reduceSubstituents(ranking)
    },
    _links {
      _makeLinks(ranking)
    },
    symmetry {symmetry},
    centerAtom {centerAtom} 
{
  // Generate the set of unique assignments possible here
  _uniqueAssignments = UniqueAssignments::uniqueAssignments(
    AssignmentType(
      symmetry,
      _reduceNeighborCharMap(_neighborCharMap),
      _links
    ),
    false // Do NOT remove trans-spanning linked groups
  );
}

/* Private members */

std::map<AtomIndexType, unsigned> CNStereocenter::_makeNeighborSymmetryPositionMap(
  const UniqueAssignments::Assignment& assignment,
  const std::map<AtomIndexType, char> neighborCharMap
) {
  std::map<AtomIndexType, unsigned> neighborSymmetryPositionMap;

  /* First get the symmetry position mapping (char -> unsigned)
   * this is e.g. map{'A' -> vector{0,2,3}, 'B' -> vector{1}}
   */
  auto charSymmetryPositionsMap = assignment.getCharMap();

  /* assign next neighbor indices using _neighborCharMap, which stores
   * neighbor's AtomIndexType -> 'A' char mapping,
   * e.g. map{4 -> 'A', 16 -> 'B', 23 -> 'A', 26 -> 'A'}
   */
  for(const auto& indexCharPair: neighborCharMap) {
    assert(
      !charSymmetryPositionsMap.at(
        indexCharPair.second // the current index's character, e.g. 'A'
      ).empty() // meaning there are symmetry positions left to assign
    );

    /* reference for better readability: the current character's symmetry
     * positions list:
     */
    std::vector<unsigned>& symmetryPositionsList = charSymmetryPositionsMap.at(
      indexCharPair.second // current character
    );

    // assign in the map
    neighborSymmetryPositionMap[
      indexCharPair.first
    ] = symmetryPositionsList.at( 
      0 // the first of the available symmetry positions for that char
    );

    // remove that first symmetry position
    symmetryPositionsList.erase(
      symmetryPositionsList.begin()
    );
  }

  return neighborSymmetryPositionMap;
}

std::map<AtomIndexType, char> CNStereocenter::_reduceSubstituents(
  const RankingInformation& ranking
) {
  /* ranking.sortedSubstituents = vector{vector{1, 4}, vector{2}, vector{3}};
   * -> reduce to {A, A, B, C}
   */

  std::map<AtomIndexType, char> indexSymbolMap;

  char letter = 'A';

  for(const auto& equalSet : ranking.sortedSubstituents) {
    for(const auto& index : equalSet) {
      indexSymbolMap[index] = letter;
    }

    ++letter;
  }

  return indexSymbolMap;
}

std::vector<char> CNStereocenter::_reduceNeighborCharMap(
  const std::map<
    AtomIndexType,
    char
  >& neighborCharMap
) const {
  std::vector<char> ligandSymbols;

  // Add every mapped char to a vector
  for(const auto& indexCharPair: neighborCharMap) {
    ligandSymbols.push_back(indexCharPair.second);
  }

  // sort it
  std::sort(
    ligandSymbols.begin(),
    ligandSymbols.end()
  );

  return ligandSymbols;
}

UniqueAssignments::Assignment::LinksSetType CNStereocenter::_makeLinks(
  const RankingInformation& ranking
  /*const std::vector<AtomIndexType>& sortedIndices,
  const std::set<
    std::pair<AtomIndexType, AtomIndexType>
  >& linkedPairsSet*/
) const {
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
void CNStereocenter::assign(const unsigned& assignment) {
  assert(assignment < _uniqueAssignments.size());

  // Store current assignment
  assignmentOption = assignment;

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndexType -> unsigned).
   */
  _neighborSymmetryPositionMap = _makeNeighborSymmetryPositionMap(
    _uniqueAssignments[assignment],
    _neighborCharMap
  );
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

  // If this stereocenter is unassigned, do nothing
  if(assignmentOption) { 
    // Compare to the current _neighborCharMap and _links
    auto newNeighborCharMap = _reduceSubstituents(newRanking);

    auto newLinks = _makeLinks(newRanking);

    if(
      _neighborCharMap != newNeighborCharMap
      || _links != newLinks
    ) {
      auto newUniqueAssignments = UniqueAssignments::uniqueAssignments(
        AssignmentType(
          symmetry,
          _reduceNeighborCharMap(newNeighborCharMap),
          newLinks
        ),
        false // do NOT remove trans-spanning linked groups
      );

      // Look for a perfect match of the neighbor symmetry position map
      auto foundIter = std::find_if(
        newUniqueAssignments.begin(),
        newUniqueAssignments.end(),
        [&](const auto& assignment) -> bool {
          return (
            _neighborSymmetryPositionMap 
            == _makeNeighborSymmetryPositionMap(
              assignment,
              newNeighborCharMap
            )
          );
        }
      );

      if(foundIter != newUniqueAssignments.end()) {
        /* Overwrite everything and set the index of this assignment as the new
         * assignment
         */
        _neighborCharMap = newNeighborCharMap;
        _links = newLinks;
        _uniqueAssignments = newUniqueAssignments;
        // _neighborSymmetryPositionMap is the same, see test in find_if

        assignmentOption = foundIter - newUniqueAssignments.begin();
      } else {
        // Overwrite but mark unassigned
        _neighborCharMap = newNeighborCharMap;
        _links = newLinks;
        _uniqueAssignments = newUniqueAssignments;

        assignmentOption = boost::none;
      }
    }
  }
}

void CNStereocenter::changeSymmetry(const Symmetry::Name& symmetryName) {
  // Set a new symmetry
  symmetry = symmetryName;

  // recalculate the number of unique Assignments
  _uniqueAssignments = UniqueAssignments::uniqueAssignments(
    AssignmentType(
      symmetryName,
      _reduceNeighborCharMap(
        _neighborCharMap
      ),
      _links
    ),
    false // do NOT remove trans-spanning ligand groups
  );

  // The Stereocenter is now unassigned
  assignmentOption = boost::none;
}

void CNStereocenter::fit(const Delib::PositionCollection& positions) {
  // Extract a list of adjacent indices from stored state
  std::vector<AtomIndexType> adjacentAtoms;
  for(const auto& mapIterPair : _neighborCharMap) {
    adjacentAtoms.push_back(mapIterPair.first);
  }

  const Symmetry::Name priorSymmetry {this->symmetry};
  const boost::optional<unsigned> priorAssignment {this->assignmentOption};

  const Symmetry::Name initialSymmetry {Symmetry::Name::Linear};
  const unsigned initialAssignment = 0;
  const unsigned initialPenalty = 100;

  Symmetry::Name bestSymmetry = initialSymmetry;
  unsigned bestAssignment = initialAssignment;
  double bestPenalty = initialPenalty;
  unsigned bestAssignmentMultiplicity = 1;

  // Cycle through all symmetries
  for(const auto& symmetryName : Symmetry::allNames) {
    // Skip any Symmetries of different size
    if(Symmetry::size(symmetryName) != Symmetry::size(symmetry)) {
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
                centerAtom,
                k
              ) - angle(
                i,
                centerAtom,
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
                  centerAtom
                ),
                DelibHelpers::getDistance( // j-k 1-2 distance from positions
                  positions,
                  centerAtom,
                  k
                ),
                angle( // idealized Stereocenter angle
                  i,
                  centerAtom,
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
    this->assignmentOption = priorAssignment;
  } else {
    // Set to best fit
    changeSymmetry(bestSymmetry);

    /* How to handle multiplicity? 
     * Current policy: If there is multiplicity, warn and do not assign
     */
    if(bestAssignmentMultiplicity > 1) {
      assignmentOption = boost::none;
    } else {
      assignmentOption = bestAssignment;
    }
  }
}

/* Information */
double CNStereocenter::angle(
  const AtomIndexType& i,
  const AtomIndexType& j, 
  const AtomIndexType& k
) const {
  assert(j == centerAtom);

  /* j is pracitcally unused here because in the Symmetry angleFunctions, the
   * middle atom is implicit, it has no symmetryPosition number. j is however
   * needed in the interface because of EZStereocenter, where it is important
   * to specify which atom is the intermediate (there are two possibilities)
   */
  return Symmetry::angleFunction(symmetry)(
    _neighborSymmetryPositionMap.at(i),
    _neighborSymmetryPositionMap.at(k)
  );
}

boost::optional<unsigned> CNStereocenter::assigned() const {
  return assignmentOption;
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
      _neighborSymmetryPositionMap
    );

    // Get list of tetrahedra from symmetry
    auto tetrahedraList = Symmetry::tetrahedra(symmetry);

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

          return centerAtom;
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
    + std::to_string(centerAtom) + " ("s + Symmetry::name(symmetry) +", "s;

  std::string charRep;
  for(const auto& iterPair : _neighborCharMap) {
    charRep += iterPair.second;
  }

  std::sort(
    charRep.begin(),
    charRep.end()
  );

  returnString += charRep + "): "s;

  if(assignmentOption) {
    returnString += std::to_string(assignmentOption.value());
  } else {
    returnString += "u";
  }

  returnString += "/"s + std::to_string(numAssignments());

  return returnString;
}

std::string CNStereocenter::rankInfo() const {
  // TODO revisit as soon as pseudo-asymmetry is added

  return (
    "CN-"s + std::to_string(static_cast<unsigned>(symmetry)) 
    + "-"s + std::to_string(numAssignments())
    + "-"s + (
      assigned() 
      ? std::to_string(assigned().value()) 
      : "u"s
    )
  );
}

std::set<AtomIndexType> CNStereocenter::involvedAtoms() const {
  return {centerAtom};
}

unsigned CNStereocenter::numAssignments() const {
  return _uniqueAssignments.size();
}

Type CNStereocenter::type() const {
  return Type::CNStereocenter;
}

bool CNStereocenter::operator == (const CNStereocenter& other) const {
  return (
    symmetry == other.symmetry
    && centerAtom == other.centerAtom
    && _neighborCharMap == other._neighborCharMap
    && _uniqueAssignments.size() == other._uniqueAssignments.size()
    && assignmentOption == other.assignmentOption
  );
}

bool CNStereocenter::operator < (const CNStereocenter& other) const {
  using TemplateMagic::componentSmaller;
  
  /* Sequentially compare individual components, comparing assignments last
   * if everything else matches
   */
  return componentSmaller(
    centerAtom,
    other.centerAtom
  ).value_or(
    componentSmaller(
      _neighborCharMap,
      other._neighborCharMap
    ).value_or(
      componentSmaller(
        _uniqueAssignments.size(),
        other._uniqueAssignments.size()
      ).value_or(
        componentSmaller(
          symmetry,
          other.symmetry
        ).value_or(
          assignmentOption < other.assignmentOption
          // NOTE: boost::none is smaller than 0 in this ordering
        )
      )
    )
  );
}

} // namespace Stereocenters

} // namespace MoleculeManip
