#include "EZStereocenter.h"

using namespace MoleculeManip;
using namespace MoleculeManip::Stereocenters;
using namespace std::string_literals;

/* Constructors */
EZStereocenter::EZStereocenter(
  const AtomIndexType& firstCenter,
  const std::vector<AtomIndexType>& firstCenterRanking,
  const AtomIndexType& secondCenter,
  const std::vector<AtomIndexType>& secondCenterRanking
) {
  assert(firstCenterRanking.size() == 2 && secondCenterRanking.size() == 2);

  // Initialize the array
  _indicesAndRank[0] = firstCenter;
  _indicesAndRank[1] = firstCenterRanking[0];
  _indicesAndRank[2] = firstCenterRanking[1];
  _indicesAndRank[3] = secondCenter;
  _indicesAndRank[4] = secondCenterRanking[0];
  _indicesAndRank[5] = secondCenterRanking[1];
}

/* Private members */
std::vector<
  std::array<AtomIndexType, 4>
> EZStereocenter::_equalPriorityDihedralSequences() const {
  return {
    { // High to High
      _indicesAndRank[1],
      _indicesAndRank[0],
      _indicesAndRank[3],
      _indicesAndRank[4]
    },{ // Low to Low
      _indicesAndRank[2],
      _indicesAndRank[0],
      _indicesAndRank[3],
      _indicesAndRank[5]
    }
  };
}

std::vector<
  std::array<AtomIndexType, 4>
> EZStereocenter::_differentPriorityDihedralSequences() const {
  return {
    { // High to Low
      _indicesAndRank[1],
      _indicesAndRank[0],
      _indicesAndRank[3],
      _indicesAndRank[5]
    },{ // Low to High
      _indicesAndRank[2],
      _indicesAndRank[0],
      _indicesAndRank[3],
      _indicesAndRank[4]
    }
  };
}

std::pair<double, double> EZStereocenter::_cisLimits() const {
  return {0, _dihedralAngleVariance};
}

std::pair<double, double> EZStereocenter::_transLimits() const {
  return {180 - _dihedralAngleVariance, 180};
}

/* Public members */
double EZStereocenter::angle(
  const AtomIndexType& i __attribute__ ((unused)),
  const AtomIndexType& j,
  const AtomIndexType& k __attribute__ ((unused))
) const {
  assert(j == _indicesAndRank[0] || j == _indicesAndRank[3]);
  /* Little optimization here -> As long as the triple i-j-k is valid, the angle
   * is always 120°
   */

  return 120;
}

void EZStereocenter::assign(const unsigned& assignment) {
  assert(assignment < 2); // so only 0 or 1
  _isEOption = static_cast<bool>(assignment);
}

boost::optional<unsigned> EZStereocenter::assigned() const {
  if(_isEOption) {
    return static_cast<unsigned>(_isEOption.value());
  } else return {};
}

unsigned EZStereocenter::assignments() const {
  return 2;
}

std::vector<Stereocenter::ChiralityConstraintPrototype> EZStereocenter::chiralityConstraints() const {
  // Three fixed ChiralityConstraints to enforce six-atom coplanarity
  
  using Prototype = Stereocenter::ChiralityConstraintPrototype;

  std::vector<Prototype> constraints {
    Prototype(
      _indicesAndRank[1],
      _indicesAndRank[2],
      _indicesAndRank[4],
      _indicesAndRank[5],
      Stereocenter::ChiralityConstraintTarget::Flat
    ),
    Prototype(
      _indicesAndRank[0],
      _indicesAndRank[1],
      _indicesAndRank[2],
      _indicesAndRank[3],
      Stereocenter::ChiralityConstraintTarget::Flat
    ),
    Prototype(
      _indicesAndRank[0],
      _indicesAndRank[3],
      _indicesAndRank[4],
      _indicesAndRank[5],
      Stereocenter::ChiralityConstraintTarget::Flat
    )
  };

  return constraints;
}

std::vector<Stereocenter::DihedralLimits> EZStereocenter::dihedralLimits() const {
  std::vector<DihedralLimits> limits;
  // EZStereocenters can impose dihedral limits 
  if(!_isEOption)  { // unspecified
    return limits;
  } else if(_isEOption.value()) {
    // In E, equal priorities are trans
    for(const auto& dihedralSequence : _equalPriorityDihedralSequences()) {
      limits.emplace_back(
        dihedralSequence,
        _transLimits()
      );
    }
    for(const auto& dihedralSequence : _differentPriorityDihedralSequences()) {
      limits.emplace_back(
        dihedralSequence,
        _cisLimits()
      );
    }

    return limits;
  } else {
    // In Z, equal priorities are cis
    for(const auto& dihedralSequence : _equalPriorityDihedralSequences()) {
      limits.emplace_back(
        dihedralSequence,
        _cisLimits()
      );
    }
    for(const auto& dihedralSequence : _differentPriorityDihedralSequences()) {
      limits.emplace_back(
        dihedralSequence,
        _transLimits()
      );
    }

    return limits;
  }
}

std::set<AtomIndexType> EZStereocenter::involvedAtoms() const {
  return {_indicesAndRank[0], _indicesAndRank[3]};
}

std::string EZStereocenter::type() const {
  std::string returnString =  "EZStereocenter on ("s 
    + std::to_string(_indicesAndRank[0]) + ", "s 
    + std::to_string(_indicesAndRank[3]) + "), "s;

  if(_isEOption) returnString += (_isEOption.value())
    ? "E"s
    : "Z"s;
  else returnString += "u";

  return returnString;
}
