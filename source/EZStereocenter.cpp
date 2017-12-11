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
) : _leftCenter {firstCenter},
    _rightCenter {secondCenter},
    _leftRanking {firstCenterRanking},
    _rightRanking {secondCenterRanking}
{
  /* Valid ranking inputs
   * - One set of one index (single substituent on double bond side)
   * - One set of two indices each (equal substituents on side)
   * - Two sets of one index each (different substituents on side)
   */
#ifndef NDEBUG
  unsigned numFirst = _numIndices(firstCenterRanking);
  unsigned numSecond = _numIndices(secondCenterRanking);

  assert(1 <= numFirst && numFirst <= 2);
  assert(1 <= numSecond && numSecond <= 2);
#endif
}

unsigned EZStereocenter::_numIndices(const RankingInformation& ranking) {
  unsigned countElements = 0;

  for(const auto& set : ranking.sortedSubstituents) {
    countElements += set.size();
  }

  return countElements;
}

/* Private members */
std::vector<
  std::array<AtomIndexType, 4>
> EZStereocenter::_equalPriorityDihedralSequences() const {
  std::vector<
    std::array<AtomIndexType, 4>
  > sequences {
    { // High to High always exists
      _leftHighPriority(),
      _leftCenter,
      _rightCenter,
      _rightHighPriority()
    }
  };

  // Low to low exists only if both exist
  if(_numIndices(_leftRanking) == 2 && _numIndices(_rightRanking) == 2) {
    sequences.emplace_back(
      std::array<AtomIndexType, 4> {
        _leftLowPriority(),
        _leftCenter,
        _rightCenter,
        _rightLowPriority()
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

  if(_numIndices(_leftRanking) == 2) {
    sequences.emplace_back(
      std::array<AtomIndexType, 4> {
        _leftLowPriority(),
        _leftCenter,
        _rightCenter,
        _rightHighPriority()
      }
    );
  }

  if(_numIndices(_rightRanking) == 2) {
    sequences.emplace_back(
      std::array<AtomIndexType, 4> {
        _leftHighPriority(),
        _leftCenter,
        _rightCenter,
        _rightLowPriority()
      }
    );
  }

  return sequences;
}

/* Public members */
/* Modification */
void EZStereocenter::addSubstituent(
  const AtomIndexType& center,
  const RankingInformation& centerRanking
) {
  /* The center must be either left or right. Any call that doesn't match
   * either is malformed
   */
  assert(center == _leftCenter || center == _rightCenter);

  if(center == _leftCenter) {
    assert(
      _numIndices(_leftRanking) == 1
      && _numIndices(centerRanking) == 2
    );

    /* In case the high priority isn't the old one anymore, we have to flip the
     * current assignment (provided there is one)
     */
    if(
      numAssignments() == 2
      && _isZOption
      && _leftHighPriority() != centerRanking.sortedSubstituents.back().back()
    ) {
      _isZOption = !_isZOption.value();
    }

    // Overwrite the remaining state
    _leftRanking = centerRanking;
  } else if(center == _rightCenter) {
    assert(
      _numIndices(_rightRanking) == 1
      && _numIndices(centerRanking) == 2
    );

    /* In case the high priority isn't the old one anymore, we have to flip the
     * current assignment (provided there is one)
     */
    if(
      numAssignments() == 2
      && _isZOption
      && _rightHighPriority() != centerRanking.sortedSubstituents.back().back()
    ) {
      _isZOption = !_isZOption.value();
    }

    // Overwrite the remaining state
    _rightRanking = centerRanking;
  } 
}

void EZStereocenter::assign(const boost::optional<unsigned>& assignment) {
  if(assignment) {
    assert(assignment.value() < numAssignments()); 
    _isZOption = static_cast<bool>(assignment.value());
  } else {
    _isZOption = boost::none;
  }
}

void EZStereocenter::fit(const Delib::PositionCollection& positions) {
  /* The only sequences that we can be sure of exist is the singular one from
   * equalPriorityDihedralSequences(). There may be two of those, though, in
   * case we have four substituents, so consider that too if it exists.
   */

  const double zPenalty = TemplateMagic::sum(
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

  const double ePenalty = TemplateMagic::sum(
    TemplateMagic::map(
      _equalPriorityDihedralSequences(),
      [&](const std::array<AtomIndexType, 4>& indices) -> double {
        return M_PI - std::fabs(
          DelibHelpers::getDihedral(
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
    _isZOption = boost::none;
  } else {
    _isZOption = zPenalty < ePenalty;
  }
}

void EZStereocenter::propagateGraphChange(
  const RankingInformation& firstCenterRanking,
  const RankingInformation& secondCenterRanking
) {
  // Assume parts of the state of EZStereocenter
  assert(numAssignments() == 2);

#ifndef NDEBUG
  unsigned numFirst = _numIndices(firstCenterRanking);
  unsigned numSecond = _numIndices(secondCenterRanking);

  assert(1 <= numFirst && numFirst <= 2);
  assert(1 <= numSecond && numSecond <= 2);
#endif

  /* As long as the class state isn't overwritten yet, decide which assignment
   * the final stereocenter will have. If the index for the high-priority index
   * doesn't match the stored one, we have to flip the state once. If this 
   * occurs at both ends, then we don't have to do anything.
   */
  bool flipAssignment = false;
  if(firstCenterRanking.sortedSubstituents.back().back() != _leftHighPriority()) {
    // Negate flipAssignment
    flipAssignment = !flipAssignment;
  }

  if(secondCenterRanking.sortedSubstituents.back().back() != _rightHighPriority()) {
    // Negate flipAssignment
    flipAssignment = !flipAssignment;
  }

  // Now we can overwrite class state
  _leftRanking = firstCenterRanking;
  _rightRanking = secondCenterRanking;

  // Change the assignment if assigned and we determined that we need to negate
  if(
    numAssignments() == 2
    && _isZOption
    && flipAssignment
  ) {
    _isZOption = !_isZOption.value();
  }
}

void EZStereocenter::propagateVertexRemoval(const AtomIndexType& removedIndex) {
  auto updateIndex = [&removedIndex](AtomIndexType& index) -> void {
    if(index > removedIndex) {
      --index;
    } else if(index == removedIndex) {
      index = std::numeric_limits<AtomIndexType>::max();
    }
  };

  auto newIndex = [&removedIndex](const AtomIndexType& index) -> AtomIndexType {
    assert(index != removedIndex);

    if(index > removedIndex) {
      return index - 1;
    }

    if(index == removedIndex) {
      return std::numeric_limits<AtomIndexType>::max();
    }

    return index;
  };

  assert(_leftCenter != removedIndex);
  assert(_rightCenter != removedIndex);

  updateIndex(_leftCenter);
  updateIndex(_rightCenter);

  auto updateRanking = [&](RankingInformation& ranking) {
    for(auto& equalPrioritySet : ranking.sortedSubstituents) {
      for(auto& index : equalPrioritySet) {
        updateIndex(index);
      }
    }

    RankingInformation::LinksType newLinks;
    for(const auto& linkPair : ranking.linkedPairs) {
      newLinks.emplace(
        newIndex(linkPair.first),
        newIndex(linkPair.second)
      );
    }
    ranking.linkedPairs = newLinks;
  };

  updateRanking(_leftRanking);
  updateRanking(_rightRanking);
}

void EZStereocenter::removeSubstituent(
  const AtomIndexType& center,
  const AtomIndexType& which
) {
  assert(center == _leftCenter || center == _rightCenter);

  auto dropWhich = [&which](RankingInformation& ranking) {
    for(auto& equalPrioritySet : ranking.sortedSubstituents) {
      TemplateMagic::inplaceRemoveIf(
        equalPrioritySet,
        [&which](const auto& index) -> bool {
          return index == which;
        }
      );
    }

    // There cannot be any linked pairs for a single substituent
    ranking.linkedPairs = {};
  };

  if(center == _leftCenter) {
    // Drop which from that ranking
    dropWhich(_leftRanking);
  } else if(center == _rightCenter) {
    dropWhich(_rightRanking);
  }
}


/* Information */
double EZStereocenter::angle(
  const AtomIndexType& i __attribute__ ((unused)),
  const AtomIndexType& j __attribute__ ((unused)),
  const AtomIndexType& k __attribute__ ((unused))
) const {
  assert(j == _leftCenter || j == _rightCenter);
  /* Little optimization here -> As long as the triple i-j-k is valid, the angle
   * is always 120°
   */

  return ConstexprMagic::Math::toRadians<double>(120);
}

boost::optional<unsigned> EZStereocenter::assigned() const {
  if(_isZOption) {
    return static_cast<unsigned>(
      _isZOption.value()
    );
  }

  return {};
}

unsigned EZStereocenter::numAssignments() const {
  /* Determine whether there can be two assignments or not. There is only one
   * assignment in the case that on either side, there are two equal
   * substituents
   */
  if(
    (
      _numIndices(_leftRanking) == 2
      && _leftRanking.sortedSubstituents.size() == 1
    ) || (
      _numIndices(_rightRanking) == 2
      && _rightRanking.sortedSubstituents.size() == 1
    )
  ) {
    return 1;
  }

  return 2;
}

std::vector<ChiralityConstraintPrototype> EZStereocenter::chiralityConstraints() const {
  // Three fixed ChiralityConstraints to enforce six-atom coplanarity
  
  std::vector<ChiralityConstraintPrototype> constraints {
    {
      std::array<AtomIndexType, 4> {
        _leftHighPriority(),
        _leftCenter,
        _rightCenter,
        _rightHighPriority(),
      },
      ChiralityConstraintTarget::Flat
    }
  };

  if(_numIndices(_leftRanking) == 2) {
    constraints.emplace_back(
      std::array<AtomIndexType, 4> {
        _leftHighPriority(),
        _leftCenter,
        _rightCenter,
        _leftLowPriority()
      },
      ChiralityConstraintTarget::Flat
    );
  }

  if(_numIndices(_rightRanking) == 2) {
    constraints.emplace_back(
      std::array<AtomIndexType, 4> {
        _rightHighPriority(),
        _rightCenter,
        _leftCenter,
        _rightLowPriority()
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
   * when _numAssignments == 1, it is irrelevant which state _isZOption is in,
   * so long as it is set (which it is, already in the constructor). Since both
   * possibilities of considering high-low on one side are equivalent, we just
   * consider the first in all cases. As long as the ranking algorithm works
   * as it should, there should be no issues.
   *
   * So we just make this dependent on the current _isZOption settings.
   */

  // EZStereocenters can impose dihedral limits 
  if(_isZOption && _isZOption.value()) {
    return _cisDihedralLimits();
  }

  if(_isZOption && !_isZOption.value()) {
    return _transDihedralLimits();
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

  if(_numIndices(_leftRanking) == 2) {
    returnString += "[H:"s + std::to_string(_leftHighPriority())
      + ", L:"s + std::to_string(_leftLowPriority()) + "], "s;
  } else {
    returnString += std::to_string(_leftHighPriority()) +", "s;
  }

  returnString +=  std::to_string(_leftCenter) + ", "s 
    + std::to_string(_rightCenter) + ", "s;

  if(_numIndices(_rightRanking) == 2) {
    returnString += "[H:"s + std::to_string(_rightHighPriority())
      + ", L:"s + std::to_string(_rightLowPriority()) + "]"s;
  } else {
    returnString += std::to_string(_rightHighPriority());
  }

  if(numAssignments() == 1) {
    returnString += ". Is non-stereogenic.";
  } else if(numAssignments() == 2) {
    returnString += ". Is "s;
    if(_isZOption) {
      returnString += (_isZOption.value())
        ? "Z"s
        : "E"s;
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

std::vector<AtomIndexType> EZStereocenter::involvedAtoms() const {
  // Return a standard form of smaller first
  return {
    std::min(_leftCenter, _rightCenter),
    std::max(_leftCenter, _rightCenter)
  };
}

Type EZStereocenter::type() const {
  return Type::EZStereocenter;
}

/* Operators */
bool EZStereocenter::operator == (const EZStereocenter& other) const {
  return (
    _leftCenter == other._leftCenter
    && _leftRanking.sortedSubstituents == other._leftRanking.sortedSubstituents
    && _rightRanking.sortedSubstituents == other._rightRanking.sortedSubstituents
    && _leftRanking.linkedPairs == other._leftRanking.linkedPairs
    && _rightRanking.linkedPairs == other._rightRanking.linkedPairs
    && _isZOption == other._isZOption
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
        _leftRanking.sortedSubstituents,
        other._leftRanking.sortedSubstituents
      ).value_or(
        componentSmaller(
          _rightRanking.sortedSubstituents,
          other._rightRanking.sortedSubstituents
        ).value_or(
          _isZOption < other._isZOption
        )
      )
    )
  );
}

// Static data
const double EZStereocenter::_dihedralAngleVariance = ConstexprMagic::Math::toRadians<double>(5);

} // namespace Stereocenters

} // namespace MoleculeManip
