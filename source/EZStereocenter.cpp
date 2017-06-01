#include "EZStereocenter.h"
#include "template_magic/Random.h"
#include "constexpr_magic/Math.h"
#include "DelibHelpers.h"

using namespace MoleculeManip;
using namespace MoleculeManip::Stereocenters;
using namespace std::string_literals;

/* Constructors */
EZStereocenter::EZStereocenter(
  const AtomIndexType& firstCenter,
  const std::vector<AtomIndexType>& firstCenterRanking,
  const IndexPairsSet& firstCenterEqualPairs,
  const AtomIndexType& secondCenter,
  const std::vector<AtomIndexType>& secondCenterRanking,
  const IndexPairsSet& secondCenterEqualPairs
) {
  assert(
    1 <= firstCenterRanking.size()
    && firstCenterRanking.size() <= 2
    && 1 <= secondCenterRanking.size() 
    && secondCenterRanking.size() <= 2
  );

  _leftCenter = firstCenter;
  _leftHighPriority = firstCenterRanking.front();
  _rightCenter = secondCenter;
  _rightHighPriority = secondCenterRanking.front();

  if(firstCenterRanking.size() == 2) {
    _leftLowPriority = firstCenterRanking.back();
  }
  if(secondCenterRanking.size() == 2) {
    _rightLowPriority = secondCenterRanking.back();
  }

  // Determine whether there can be two assignments or not
  if(firstCenterEqualPairs.empty() && secondCenterEqualPairs.empty()) {
    _numAssignments = 2;
  } else {
    _numAssignments = 1;
    _isEOption = false;
  }
}

/* Private members */
std::vector<
  std::array<AtomIndexType, 4>
> EZStereocenter::_equalPriorityDihedralSequences() const {
  std::vector<
    std::array<AtomIndexType, 4>
  > sequences {
    { // High to High always exists
      _leftHighPriority,
      _leftCenter,
      _rightCenter,
      _rightHighPriority
    }
  };

  // Low to low exists only if both exist
  if(_leftLowPriority && _rightLowPriority) {
    sequences.emplace_back(
      std::array<AtomIndexType, 4> {
        _leftLowPriority.value(),
        _leftCenter,
        _rightCenter,
        _rightLowPriority.value()
      }
    );
  }

  return sequences;
}

std::vector<
  std::array<AtomIndexType, 4>
> EZStereocenter::_differentPriorityDihedralSequences() const {
  
  std::vector<
    std::array<AtomIndexType, 4>
  > sequences;

  if(_leftLowPriority) {
    sequences.emplace_back(
      std::array<AtomIndexType, 4> {
        _leftLowPriority.value(),
        _leftCenter,
        _rightCenter,
        _rightHighPriority
      }
    );
  }

  if(_rightLowPriority) {
    sequences.emplace_back(
      std::array<AtomIndexType, 4> {
        _leftHighPriority,
        _leftCenter,
        _rightCenter,
        _rightLowPriority.value()
      }
    );
  }

  return sequences;
}

/* Public members */
/* Modification */
void EZStereocenter::assign(const unsigned& assignment) {
  assert(assignment < _numAssignments); 

  _isEOption = static_cast<bool>(assignment);
}

void EZStereocenter::fit(const Delib::PositionCollection& positions) {
  /* The only sequences that we can be sure of exist is the singular one from
   * equalPriorityDihedralSequences(). There may be two of those, though, in
   * case we have four substituents, so consider that too if it exists.
   */

  double zPenalty = TemplateMagic::numeric::sum(
    TemplateMagic::map(
      _equalPriorityDihedralSequences(),
      [&](const std::array<AtomIndexType, 4>& indices) -> double {
        return std::fabs(
          DelibHelpers::getDihedral(
            positions,
            indices
          )
        );
      }
    )
  );

  double ePenalty = TemplateMagic::numeric::sum(
    TemplateMagic::map(
      _equalPriorityDihedralSequences(),
      [&](const std::array<AtomIndexType, 4>& indices) -> double {
        return std::fabs(
          M_PI - DelibHelpers::getDihedral(
            positions,
            indices
          )
        );
      }
    )
  );

  if(zPenalty == ePenalty) {
    /* This means our substituents are at right angles to one another for some
     * reason or another
     */
    _isEOption = boost::none;
  } else {
    _isEOption = zPenalty > ePenalty;
  }
}

/* Information */
double EZStereocenter::angle(
  const AtomIndexType& i __attribute__ ((unused)),
  const AtomIndexType& j,
  const AtomIndexType& k __attribute__ ((unused))
) const {
  assert(j == _leftCenter || j == _rightCenter);
  /* Little optimization here -> As long as the triple i-j-k is valid, the angle
   * is always 120°
   */

  return ConstexprMagic::Math::toRadians(120);
}

boost::optional<unsigned> EZStereocenter::assigned() const {
  if(_isEOption) {
    return static_cast<unsigned>(
      _isEOption.value()
    );
  }

  return {};
}

unsigned EZStereocenter::numAssignments() const {
  return _numAssignments;
}

std::vector<ChiralityConstraintPrototype> EZStereocenter::chiralityConstraints() const {
  // Three fixed ChiralityConstraints to enforce six-atom coplanarity
  
  std::vector<ChiralityConstraintPrototype> constraints {
    {
      std::array<AtomIndexType, 4> {
        _leftHighPriority,
        _leftCenter,
        _rightCenter,
        _rightHighPriority
      },
      ChiralityConstraintTarget::Flat
    }
  };

  return constraints;
}

std::vector<DihedralLimits> EZStereocenter::_cisDihedralLimits() const {
  std::vector<DihedralLimits> limits;

  // In Z, equal priorities are cis
  for(const auto& dihedralSequence : _equalPriorityDihedralSequences()) {
    limits.emplace_back(
      dihedralSequence,
      std::pair<double, double> {0, _dihedralAngleVariance} // cis limit
    );
  }

  for(const auto& dihedralSequence : _differentPriorityDihedralSequences()) {
    limits.emplace_back(
      dihedralSequence,
      std::pair<double, double> {M_PI - _dihedralAngleVariance, M_PI} // trans
    );
  }

  return limits;
}

std::vector<DihedralLimits> EZStereocenter::_transDihedralLimits() const {
  std::vector<DihedralLimits> limits;

  // In E, equal priorities are trans
  for(const auto& dihedralSequence : _equalPriorityDihedralSequences()) {
    limits.emplace_back(
      dihedralSequence,
      std::pair<double, double> {M_PI - _dihedralAngleVariance, M_PI} // trans
    );
  }

  for(const auto& dihedralSequence : _differentPriorityDihedralSequences()) {
    limits.emplace_back(
      dihedralSequence,
      std::pair<double, double> {0, _dihedralAngleVariance} // cis limit
    );
  }

  return limits;
}

std::vector<DihedralLimits> EZStereocenter::dihedralLimits() const {
  /* Whether we have _numAssignments == 1 or 2 can be completely ignored because
   * when _numAssignments == 1, it is irrelevant which state _isEOption is in,
   * so long as it is set (which it is, already in the constructor). Since both
   * possibilities of considering high-low on one side are equivalent, we just
   * consider the first in all cases. As long as the ranking algorithm works
   * as it should, there should be no issues.
   *
   * So we just make this dependent on the current _isEOption settings:
   */

  // EZStereocenters can impose dihedral limits 
  if(_isEOption && _isEOption.value()) {
    return _transDihedralLimits();
  }

  if(_isEOption && !_isEOption.value()) {
    return _cisDihedralLimits();
  }

  return {};
}

std::string EZStereocenter::info() const {
  std::string returnString =  "EZStereocenter on ("s 
    + std::to_string(_leftCenter) + ", "s 
    + std::to_string(_rightCenter) + "), "s;

  if(_isEOption) {
    returnString += (_isEOption.value())
      ? "E"s
      : "Z"s;
  } else {
    returnString += "u";
  }

  return returnString;
}

std::set<AtomIndexType> EZStereocenter::involvedAtoms() const {
  return {_leftCenter, _rightCenter};
}

Type EZStereocenter::type() const {
  return Type::EZStereocenter;
}

/* Operators */
bool EZStereocenter::operator == (const EZStereocenter& other) const {
  return (
    _leftCenter == other._leftCenter
    && _leftHighPriority == other._leftHighPriority
    && _rightCenter == other._rightCenter
    && _rightHighPriority == other._rightHighPriority
    && _leftLowPriority == other._leftLowPriority
    && _numAssignments == other._numAssignments
    && _isEOption == other._isEOption
  );
}

// Static data
const double EZStereocenter::_dihedralAngleVariance = ConstexprMagic::Math::toRadians(5);
