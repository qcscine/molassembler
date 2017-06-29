#include "CNStereocenter.h"
#include "geometry_assignment/GenerateUniques.h"
#include "template_magic/TemplateMagic.h"
#include "StdlibTypeAlgorithms.h"
#include "DelibHelpers.h"
#include "CommonTrig.h"
#include "Log.h"

#include <iomanip>

namespace MoleculeManip {

namespace Stereocenters {

/* Constructors */
CNStereocenter::CNStereocenter(
  // The symmetry of this Stereocenter
  const Symmetry::Name& symmetry,
  // The atom this Stereocenter is centered on
  const AtomIndexType& centerAtom,
  // A partially ordered list of substituents, low to high 
  const std::vector<AtomIndexType>& partiallySortedSubstituents,
  // A set of pairs denoting which substituents are equal priority
  const std::set<
    std::pair<AtomIndexType, AtomIndexType>
  >& equalPairsSet
) : symmetry(symmetry),
    centerAtom(centerAtom) {

  // Reduce the substituents to a character map
  _neighborCharMap = _reduceSubstituents(
    partiallySortedSubstituents,
    equalPairsSet
  );

  // Generate the set of unique assignments possible here
  _uniqueAssignments = UniqueAssignments::uniqueAssignments(
    AssignmentType(
      symmetry,
      _reduceNeighborCharMap(
        _neighborCharMap
      )
    )
  );

  /* auto reducedMap = _reduceNeighborCharMap(_neighborCharMap);
  std::cout << "Reduced char map: vec{";
  for(const auto& character : reducedMap) {
    std::cout << character;
    if(character != reducedMap.back()) std::cout << ", ";
  }
  std::cout << "} in symmetry " << Symmetry::name(symmetry) 
    << " -> " << _uniqueAssignments.size() << " assignments" << std::endl;*/
}

/* Private members */
std::map<
  AtomIndexType,
  char
> CNStereocenter::_reduceSubstituents(
  const std::vector<AtomIndexType>& rankedSubstituentNextAtomIndices,
  const std::set<
    std::pair<
      AtomIndexType,
      AtomIndexType
    >
  >& equalSubstituentPairsSet
) const {
  /* the algorithm returns pairs of ligands that are equal (due to some 
   * custom sorting function shenanigans, which is binary). We want to 
   * condense that information into sets of equal ligands, so we restructure
   * the set of pairs to a vector of non-overlapping sets:
   * e.g. set{pair{1, 3}, pair{1, 4}, pair{2, 5}} 
   *  -> vector{set{1, 3, 4}, set{2, 5}}
   */
  auto setsVector = StdlibTypeAlgorithms::makeIndividualSets(
    equalSubstituentPairsSet
  );

  // Add lone substituents to setsVector
  for(const auto& index: rankedSubstituentNextAtomIndices) {
    // if the current substituent index is not in any of the sets
    if(!std::accumulate(
      setsVector.begin(),
      setsVector.end(),
      false,
      [&index](const bool& carry, const std::set<AtomIndexType>& set) {
        return (
          carry
          || set.count(index) == 1
        );
      }
    )) {
      // add a single-atom set
      setsVector.emplace_back(
        std::set<AtomIndexType>{index}
      );
    }
  }

  /* so now we have e.g.
   * setsVector = vector{set{1, 4}, set{2}, set{3}};
   * rankedSubstituentNextAtomIndices = vector{ 2, 1, 4, 3};
   * -> reduce to {A, B, B, C}
   */

  // create a mapping between indices and ligand symbols
  std::map<
    AtomIndexType,
    char
  > indexSymbolMap;

  const char initialChar = 'A';
  for(const auto& index : rankedSubstituentNextAtomIndices) {
    // find position in setsVector

    auto findIter = std::find_if(
      setsVector.begin(),
      setsVector.end(),
      [&index](const auto& set) -> bool {
        return set.count(index) > 0;
      }
    );

    assert(findIter != setsVector.end()); // would indicate an error above

    unsigned posInSetsVector = findIter - setsVector.begin();

    // Take advantage of implicit type conversions:
    //               char =        char +        unsigned
    indexSymbolMap[index] = initialChar + posInSetsVector;
  }


  return indexSymbolMap;

  /* TODO no use of connectivity information as of yet to determine whether 
   * ligands are bridged!
   */
}

std::vector<char> CNStereocenter::_reduceNeighborCharMap(
  const std::map<
    AtomIndexType,
    char
  >& neighborCharMap
) {
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

/* Modification */
void CNStereocenter::assign(const unsigned& assignment) {
  assert(assignment < _uniqueAssignments.size());

  // Store current assignment
  assignmentOption = assignment;

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndexType -> unsigned).
   *
   * First get the symmetry position mapping (char -> unsigned)
   * this is e.g. map{'A' -> vector{0,2,3}, 'B' -> vector{1}}
   */
  auto charSymmetryPositionsMap = _uniqueAssignments[
    assignment
  ].getCharMap();

  /* assign next neighbor indices using _neighborCharMap, which stores
   * neighbor's AtomIndexType -> 'A' char mapping,
   * e.g. map{4 -> 'A', 16 -> 'B', 23 -> 'A', 26 -> 'A'}
   */
  for(const auto& indexCharPair: _neighborCharMap) {
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
    _neighborSymmetryPositionMap[
      indexCharPair.first
    ] = symmetryPositionsList.at( 
      0 // the first of the available symmetry positions for that char
    );

    // remove that first symmetry position
    symmetryPositionsList.erase(
      symmetryPositionsList.begin()
    );
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
      )
    )
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

  Symmetry::Name bestSymmetry = Symmetry::Name::Linear;
  unsigned bestAssignment = 0;
  double bestPenalty = 100;
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

      const double angleDeviations = TemplateMagic::numeric::sum(
        TemplateMagic::allPairsMap(
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

      const double oneThreeDistanceDeviations = TemplateMagic::numeric::sum(
        TemplateMagic::allPairsMap(
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
        : TemplateMagic::numeric::sum(
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
  /* This is a bit weird, I thought the boost::optional<T>::operator ==
   * emcompassed all this. Maybe the behavior is different in a more recent
   * Boost version, or just switching to std::optional with C++17 will allow
   * straight-up use of operator ==.
   */
  bool assignmentsSame = (
    (
      assignmentOption 
      && other.assignmentOption 
      && assignmentOption.value() == other.assignmentOption.value()
    ) || (!assignmentOption && !other.assignmentOption)
  );

  /*if(!assignmentsSame) {
    std::cout << "Different assignment" << std::endl;
    std::cout << "This:" << assignment.value_or(50) << ", Other:" 
      << other.assignment.value_or(50) << std::endl;
    std::cout << std::boolalpha << "This is assigned: " 
      << static_cast<bool>(assignment) << ", Other is assigned: " 
      << static_cast<bool>(other.assignment) << std::endl;
  }
  if(symmetry != other.symmetry) std::cout << "Symmetry is different." << std::endl;
  if(centerAtom != other.centerAtom) std::cout << "Central atom is different." << std::endl;
  if(_neighborCharMap != other._neighborCharMap) std::cout << "Char map is different." << std::endl;
  if(_uniqueAssignments.size() != other._uniqueAssignments.size()) std::cout << "Different #assignments." << std::endl;*/

  return (
    symmetry == other.symmetry
    && centerAtom == other.centerAtom
    && _neighborCharMap == other._neighborCharMap
    && _uniqueAssignments.size() == other._uniqueAssignments.size()
    && assignmentsSame
  );
}

} // namespace Stereocenters

} // namespace MoleculeManip
