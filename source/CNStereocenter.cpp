#include "CNStereocenter.h"
#include "steric_uniqueness/GenerateUniques.h"
#include "template_magic/templateMagic.h"

namespace MoleculeManip {

namespace Stereocenters {

/* Constructors */
CNStereocenter::CNStereocenter(
  // The symmetry of this Stereocenter
  const Symmetry::Name& symmetry,
  // The atom this Stereocenter is centered on
  const AtomIndexType& centerAtom,
  // A partially ordered list of substituents, low to high 
  const std::vector<AtomIndexType> partiallySortedSubstituents,
  // A set of pairs denoting which substituents are equal priority
  const std::set<
    std::pair<AtomIndexType, AtomIndexType>
  > equalPairsSet
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
      setsVector.push_back(
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
void CNStereocenter::assign(const unsigned& passAssignment) {
  assert(passAssignment < _uniqueAssignments.size());

  // Store current assignment
  assignment = passAssignment;

  /* save a mapping of next neighbor indices to symmetry positions after
   * assigning (AtomIndexType -> unsigned).
   *
   * First get the symmetry position mapping (char -> unsigned)
   * this is e.g. map{'A' -> vector{0,2,3}, 'B' -> vector{1}}
   */
  auto charSymmetryPositionsMap = _uniqueAssignments[
    passAssignment
  ].getCharMap();

  /* assign next neighbor indices using _neighborCharMap, which stores
   * neighbor's AtomIndexType -> 'A' char mapping,
   * e.g. map{4 -> 'A', 16 -> 'B', 23 -> 'A', 26 -> 'A'}
   */
  for(const auto& indexCharPair: _neighborCharMap) {
    assert(
      charSymmetryPositionsMap.at(
        indexCharPair.second // the current index's character, e.g. 'A'
      ).size() > 0 // meaning there are symmetry positions left to assign
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
  assignment = boost::none;
}

/* Information */
double CNStereocenter::angle(
  const AtomIndexType& i,
  const AtomIndexType& j __attribute__((unused)), 
  const AtomIndexType& k
) const {
  /* j is unused here because in the Symmetry angleFunctions, the middle atom is
   * implicit, it has no symmetryPosition number. j is however needed in the
   * interface because of EZStereocenter, where it is important to specify which
   * atom is the intermediate (there are two possibilities)
   */
  return Symmetry::angleFunction(symmetry)(
    _neighborSymmetryPositionMap.at(i),
    _neighborSymmetryPositionMap.at(k)
  );
}

boost::optional<unsigned> CNStereocenter::assigned() const {
  return assignment;
}

// TODO rename to numAssignments
unsigned CNStereocenter::numAssignments() const {
  return _uniqueAssignments.size();
}

std::vector<
  Stereocenter::ChiralityConstraintPrototype
> CNStereocenter::chiralityConstraints() const {
  std::vector<ChiralityConstraintPrototype> prototypes;
  
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
        } else {
          return centerAtom;
        }
      }
    );

    // Make a prototype from it
    prototypes.emplace_back(
      replaced,
      ChiralityConstraintTarget::Positive
    );
  }

  return prototypes;
}

std::vector<Stereocenter::DihedralLimits> CNStereocenter::dihedralLimits() const {
  return {};
}

std::string CNStereocenter::info() const {
  std::string returnString = "CNStereocenter on "s 
    + std::to_string(centerAtom) + " ("s + Symmetry::name(symmetry) +"): "s;

  if(assignment) returnString += std::to_string(assignment.value());
  else returnString += "u";

  returnString += "/"s + std::to_string(numAssignments());

  return returnString;
}

std::set<AtomIndexType> CNStereocenter::involvedAtoms() const {
  return {centerAtom};
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
    && assignment == other.assignment
  );
}

}

}
