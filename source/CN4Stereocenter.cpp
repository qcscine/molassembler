#include "CN4Stereocenter.h"

#include <math.h>
#include <Eigen/LU>
#include "CommonTrig.h"
#include "StdlibTypeAlgorithms.h"
#include "UniqueAssignments/GenerateUniques.h"

namespace MoleculeManip {

namespace Stereocenters {

std::map<
  AtomIndexType,
  char
> CN4Stereocenter::_reduceSubstituents() const {
  // OLD ––––––––––––––––––
  std::vector<AtomIndexType> rankedSubstituentNextAtomIndices;
  std::set<
    std::pair<
      AtomIndexType,
      AtomIndexType
    >
  > equalSubstituentPairsSet;

  std::tie(
    rankedSubstituentNextAtomIndices,
    equalSubstituentPairsSet
  ) = _molPtr -> rankPriority(
    _centerAtom,
    {} // nothing to exclude
  );

  /* C++17 structured binding improvement assuming tuples are still
   * convertible to pairs
  auto [
    rankedSubstituentNextAtomIndices,
    equalSubstituentPairsSet
  ] = _molPtr -> rankPriority(
    _centerAtom,
    {}
  );
  */

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
    unsigned posInSetsVector = 0;
    while(
      setsVector[posInSetsVector].count(index) == 0 
      && posInSetsVector < setsVector.size()
    ) {
      posInSetsVector++;
    }
    indexSymbolMap[index] = initialChar + posInSetsVector;
  }


  return indexSymbolMap;

  /* TODO no use of connectivity information as of yet to determine whether 
   * ligands are bridged!
   */
}

std::vector<char> CN4Stereocenter::_reduceNeighborCharMap(
  const std::map<
    AtomIndexType,
    char
  >& neighborCharMap
) {
  std::vector<char> ligandSymbols;
  for(const auto& indexCharPair: neighborCharMap) {
    ligandSymbols.push_back(indexCharPair.second);
  }

  std::sort(
    ligandSymbols.begin(),
    ligandSymbols.end()
  );

  return ligandSymbols;
}

CN4Stereocenter::CN4Stereocenter(
  const Molecule* molPtr,
  const AtomIndexType& center
) : 
  _molPtr(molPtr),
  _centerAtom (center)
{
  // save the symbolic ligand case
  _neighborCharMap = _reduceSubstituents();

  // save a list of all unique assignments for the current case
  /* TODO -> Cache these, no need to have this duplicated for every
   * stereocenter with the same ligand case
   */

  _uniqueAssignments = UniqueAssignments::uniqueAssignments(
    AssignmentType(
      _reduceNeighborCharMap(
        _neighborCharMap
      )
    )
  );

  /* save a mapping of the next neighbor's index to the position in the
   * stereocenter's symmetry representation
   */
  /* TODO assumes unspecified initially, no initialization from 3D yet
   * otherwise should save current mapping!
   */
}

void CN4Stereocenter::assign(const unsigned& assignment) {
  assert(assignment < _uniqueAssignments.size());
  if(!_assignment) { // unassigned previously
    // assign as normal
    _assignment = assignment;

    /* save a mapping of next neighbor indices to symmetry positions after
     * assigning (AtomIndexType -> uint8_t).
     *
     * First get the symmetry position mapping (char -> uint8_t)
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
        charSymmetryPositionsMap.at(
          indexCharPair.second // the current index's character, e.g. 'A'
        ).size() > 0 // meaning there are symmetry positions left to assign
      );

      /* reference for better readability: the current character's symmetry
       * positions list:
       */
      std::vector<uint8_t>& symmetryPositionsList = charSymmetryPositionsMap.at(
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

    // DONE, now _neighborSymmetryPositionMap has a mapping
  } else {
    // just assign, explicitly do NOT change the mapping to symmetry positions
    _assignment = assignment;
  }
}

std::pair<
  std::vector<DistanceConstraint>,
  std::vector<ChiralityConstraint>
> CN4Stereocenter::collectConstraints() const {

  // create the required vectors
  std::vector<DistanceConstraint> distanceConstraints;
  std::vector<ChiralityConstraint> chiralityConstraints; 

  const auto neighbors = _molPtr -> getBondedAtomIndices(_centerAtom);

  // 1-2 constraints
  for(unsigned i = 0; i < neighbors.size(); i++) {
    double a = Bond::calculateBondDistance(
      _molPtr -> getElementType(neighbors[i]),
      _molPtr -> getElementType(_centerAtom),
      _molPtr -> getBondType(
        neighbors[i],
        _centerAtom
      )
    );
    distanceConstraints.emplace_back(
      std::min(_centerAtom, neighbors[i]),
      std::max(_centerAtom, neighbors[i]), 
      a * 1.01,
      a * 0.99
    );
  }

  // Matrix for tetrahedron volume calculation later
  Eigen::Matrix<double, 5, 5> cayleyMenger;
  cayleyMenger.setZero();

  // 1-3 constraints, store distance for volume calculation later
  for(unsigned i = 0; i < neighbors.size(); i++) {
    for(unsigned j = i + 1; j < neighbors.size(); j++) {
      double a = Bond::calculateBondDistance(
        _molPtr -> getElementType(neighbors[i]),
        _molPtr -> getElementType(_centerAtom),
        _molPtr -> getBondType(
          neighbors[i],
          _centerAtom
        )
      );
      double b = Bond::calculateBondDistance(
        _molPtr -> getElementType(neighbors[j]),
        _molPtr -> getElementType(_centerAtom),
        _molPtr -> getBondType(
          neighbors[j],
          _centerAtom
        )
      );

      auto ijAngle = StdlibTypeAlgorithms::minMaxAdaptor( 
        StdlibTypeAlgorithms::makeFunction(
          PermSymmetry::Tetrahedral<UniqueAssignments::AssignmentColumn>::angle
        ),
        _neighborSymmetryPositionMap.at(neighbors[i]),
        _neighborSymmetryPositionMap.at(neighbors[j])
      );

      distanceConstraints.emplace_back(
        std::min(neighbors[i], neighbors[j]),
        std::max(neighbors[i], neighbors[j]),
        CommonTrig::lawOfCosines(
          a * 0.99,
          b * 0.99,
          ijAngle
        ),
        CommonTrig::lawOfCosines(
          a * 1.01,
          b * 1.01,
          ijAngle
        )
      );

      cayleyMenger(i + 1, j + 1) = CommonTrig::lawOfCosines(
        a,
        b,
        ijAngle
      );
    }
  }

  // top row of cayleyMenger matrix
  for(unsigned i = 1; i < 5; i++) {
    cayleyMenger(0, i) = 1;
  }

  // get the determinant
  auto determinant = static_cast<
    Eigen::Matrix<double, 5, 5>
  >(
    std::move(
      cayleyMenger.selfadjointView<Eigen::Upper>()
    )
  ).determinant();

  auto chiralityTarget = sqrt(
    determinant / 8.0
  );

  // switch target value depending on current assignment
  // TODO as soon as a full flow from reading a stereocenter to writing one
  // is completed, confirm or alter this if (==0 or ==1), one of both is 
  // correct.
  if(_assignment.value() % 2 == 0) chiralityTarget *= -1.0;

  chiralityConstraints.emplace_back(
    neighbors[0],
    neighbors[1],
    neighbors[2],
    neighbors[3],
    chiralityTarget
  );

  return make_pair(
    std::move(distanceConstraints),
    std::move(chiralityConstraints)
  );
}

} // eo namespace Stereocenters

} // eo namespace MoleculeManip
