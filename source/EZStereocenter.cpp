#include "constexpr_magic/Math.h"
#include "template_magic/Boost.h"
#include "template_magic/Containers.h"
#include "template_magic/Numeric.h"
#include "template_magic/Random.h"

#include "DelibHelpers.h"
#include "EZStereocenter.h"

namespace MoleculeManip {

namespace Stereocenters {

/* Constructors */
EZStereocenter::EZStereocenter(
  const AtomIndexType& firstCenter,
  const RankingInformation& firstCenterRanking,
  const AtomIndexType& secondCenter,
  const RankingInformation& secondCenterRanking
) {
  /* Valid ranking inputs
   * - One set of one index (single substituent on double bond side)
   * - One set of two indices each (equal substituents on side)
   * - Two sets of one index each (different substituents on side)
   */
  auto numIndices = [](const RankingInformation& ranking) -> unsigned {
    unsigned countElements = 0;

    for(const auto& set : ranking.sortedSubstituents) {
      countElements += set.size();
    }

    return countElements;
  };

  auto numFirst = numIndices(firstCenterRanking);
  auto numSecond = numIndices(secondCenterRanking);

  assert(1 <= numFirst && numFirst <= 2);
  assert(1 <= numSecond && numSecond <= 2);

  /* NOTE: high priority is assigned from the back of the ranking since the list
   * is ordered ascending
   *
   * NOTE: back-back always works out to the right selection for all valid
   * inputs.
   */
  _leftCenter = firstCenter;
  _leftHighPriority = firstCenterRanking.sortedSubstituents.back().back();
  _rightCenter = secondCenter;
  _rightHighPriority = secondCenterRanking.sortedSubstituents.back().back();

  if(numFirst == 2) {
    _leftLowPriority = firstCenterRanking.sortedSubstituents.front().front();
  }
  if(numSecond == 2) {
    _rightLowPriority = secondCenterRanking.sortedSubstituents.front().front();
  }

  // Determine whether there can be two assignments or not
  if(
    firstCenterRanking.sortedSubstituents.size() == 2
    && secondCenterRanking.sortedSubstituents.size() == 2
  ) {
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
void EZStereocenter::adaptToRankingChange(const RankingInformation& newRanking) {
  // TODO implement!
}

void EZStereocenter::assign(const unsigned& assignment) {
  assert(assignment < _numAssignments); 

  _isEOption = static_cast<bool>(assignment);
}

void EZStereocenter::fit(const Delib::PositionCollection& positions) {
  /* The only sequences that we can be sure of exist is the singular one from
   * equalPriorityDihedralSequences(). There may be two of those, though, in
   * case we have four substituents, so consider that too if it exists.
   */

  double zPenalty = TemplateMagic::sum(
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

  double ePenalty = TemplateMagic::sum(
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

  return ConstexprMagic::Math::toRadians<double>(120);
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

  if(_leftLowPriority) {
    constraints.emplace_back(
      std::array<AtomIndexType, 4> {
        _leftHighPriority,
        _leftCenter,
        _rightCenter,
        _leftLowPriority.value()
      },
      ChiralityConstraintTarget::Flat
    );
  }

  if(_rightLowPriority) {
    constraints.emplace_back(
      std::array<AtomIndexType, 4> {
        _rightHighPriority,
        _rightCenter,
        _leftCenter,
        _rightLowPriority.value()
      },
      ChiralityConstraintTarget::Flat
    );
  }

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
   * So we just make this dependent on the current _isEOption settings.
   */

  // EZStereocenters can impose dihedral limits 
  if(_isEOption && _isEOption.value()) {
    return _transDihedralLimits();
  }

  if(_isEOption && !_isEOption.value()) {
    return _cisDihedralLimits();
  }

  /* It could occur to you that the limits are defined only on the positive
   * interval [0, π], and so a corresponding trans dihedral lies on an interval
   * with some tolerance t: [t, π]. What about the negative dihedral interval?
   *
   * The way this data is used in distance geometry is merely a distance 
   * consideration depending on all bond lengths and angles involved, which is 
   * symmetric to the negative interval. By only considering the positive
   * interval, the correct distances are determined for the negative interval
   * as well.
   */

  return {};
}

std::string EZStereocenter::info() const {
  using namespace std::string_literals;

  std::string returnString =  "EZ indices "s;

  if(_leftLowPriority) {
    returnString += "[H:"s + std::to_string(_leftHighPriority)
      + ", L:"s + std::to_string(_leftLowPriority.value()) + "], "s;
  } else {
    returnString += std::to_string(_leftHighPriority) +", "s;
  }

  returnString +=  std::to_string(_leftCenter) + ", "s 
    + std::to_string(_rightCenter) + ", "s;

  if(_rightLowPriority) {
    returnString += "[H:"s + std::to_string(_rightHighPriority)
      + ", L:"s + std::to_string(_rightLowPriority.value()) + "]"s;
  } else {
    returnString += std::to_string(_rightHighPriority);
  }

  if(numAssignments() == 1) {
    returnString += ". Is non-stereogenic.";
  } else if(numAssignments() == 2) {
    returnString += ". Is "s;
    if(_isEOption) {
      returnString += (_isEOption.value())
        ? "E"s
        : "Z"s;
    } else {
      returnString += "u";
    }
  }

  return returnString;
}

std::string EZStereocenter::rankInfo() const {
  // TODO revisit as soon as pseudo-asymmetry is introduced
  using namespace std::string_literals;

  return (
    "EZ-"s + std::to_string(numAssignments())
    + "-"s + (
      assigned() 
      ? std::to_string(assigned().value()) 
      : "u"s
    )
  );
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
    && _rightLowPriority == other._rightLowPriority
    && _numAssignments == other._numAssignments
    && _isEOption == other._isEOption
  );
}

bool EZStereocenter::operator < (const EZStereocenter& other) const {
  using TemplateMagic::componentSmaller;

  return componentSmaller(
    _leftCenter,
    other._leftCenter
  ).value_or(
    componentSmaller(
    _rightCenter,
    other._rightCenter
    ).value_or(
      componentSmaller(
        _leftHighPriority,
        other._leftHighPriority
      ).value_or(
        componentSmaller(
          _rightHighPriority,
          other._rightHighPriority
        ).value_or(
          componentSmaller(
            _leftLowPriority,
            other._leftLowPriority
          ).value_or(
            componentSmaller(
              _rightLowPriority,
              other._rightLowPriority
            ).value_or(
              _isEOption < other._isEOption
            )
          )
        )
      )
    )
  );
}

// Static data
const double EZStereocenter::_dihedralAngleVariance = ConstexprMagic::Math::toRadians<double>(5);

} // namespace Stereocenters

} // namespace MoleculeManip
