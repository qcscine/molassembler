#include "temple/constexpr/ConsecutiveCompare.h"
#include "temple/constexpr/Math.h"
#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Random.h"

#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"

#include "detail/DelibHelpers.h"
#include "EZStereocenter.h"
#include "BondDistance.h"
#include "Molecule.h"

namespace molassembler {

namespace Stereocenters {

constexpr double EZStereocenter::chiralityConstraintTolerance;

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

  for(const auto& set : ranking.ligands) {
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
  RankingInformation centerRanking
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
      numStereopermutations() == 2
      && _isZOption
      && _leftHighPriority() != highPriority(centerRanking)
    ) {
      _isZOption = !_isZOption.value();
    }

    // Overwrite the remaining state
    _leftRanking = std::move(centerRanking);
  } else if(center == _rightCenter) {
    assert(
      _numIndices(_rightRanking) == 1
      && _numIndices(centerRanking) == 2
    );

    /* In case the high priority isn't the old one anymore, we have to flip the
     * current assignment (provided there is one)
     */
    if(
      numStereopermutations() == 2
      && _isZOption
      && _rightHighPriority() != highPriority(centerRanking)
    ) {
      _isZOption = !_isZOption.value();
    }

    // Overwrite the remaining state
    _rightRanking = std::move(centerRanking);
  }
}

void EZStereocenter::assign(boost::optional<unsigned> assignment) {
  if(assignment) {
    assert(assignment.value() < numStereopermutations());
    _isZOption = static_cast<bool>(assignment.value());
  } else {
    _isZOption = boost::none;
  }
}

void EZStereocenter::assignRandom() {
  assert(_isZOption == boost::none);

  _isZOption = temple::random.getSingle<bool>();
}

void EZStereocenter::fit(const AngstromWrapper& angstromWrapper) {
  /* The only sequences that we can be sure of exist is the singular one from
   * equalPriorityDihedralSequences(). There may be two of those, though, in
   * case we have four substituents, so consider that too if it exists.
   */

  const double zPenalty = temple::sum(
    temple::map(
      _equalPriorityDihedralSequences(),
      [&](const std::array<AtomIndexType, 4>& indices) -> double {
        return std::fabs(
          DelibHelpers::getDihedral(angstromWrapper.positions, indices)
        );
      }
    )
  ) + temple::sum(
    temple::map(
      _differentPriorityDihedralSequences(),
      [&](const std::array<AtomIndexType, 4>& indices) -> double {
        return M_PI - std::fabs(
          DelibHelpers::getDihedral(angstromWrapper.positions, indices)
        );
      }
    )
  );

  const double ePenalty = temple::sum(
    temple::map(
      _equalPriorityDihedralSequences(),
      [&](const std::array<AtomIndexType, 4>& indices) -> double {
        return M_PI - std::fabs(
          DelibHelpers::getDihedral(angstromWrapper.positions, indices)
        );
      }
    )
  ) + temple::sum(
    temple::map(
      _differentPriorityDihedralSequences(),
      [&](const std::array<AtomIndexType, 4>& indices) -> double {
        return std::fabs(
          DelibHelpers::getDihedral(angstromWrapper.positions, indices)
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
  RankingInformation firstCenterRanking,
  RankingInformation secondCenterRanking
) {
  // Assume parts of the state of EZStereocenter
  assert(numStereopermutations() == 2);

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
  bool flipStereopermutation = false;
  if(highPriority(firstCenterRanking) != _leftHighPriority()) {
    // Negate flipStereopermutation
    flipStereopermutation = !flipStereopermutation;
  }

  if(highPriority(secondCenterRanking) != _rightHighPriority()) {
    // Negate flipStereopermutation
    flipStereopermutation = !flipStereopermutation;
  }

  // Now we can overwrite class state
  _leftRanking = std::move(firstCenterRanking);
  _rightRanking = std::move(secondCenterRanking);

  // Change the assignment if assigned and we determined that we need to negate
  if(
    numStereopermutations() == 2
    && _isZOption
    && flipStereopermutation
  ) {
    _isZOption = !_isZOption.value();
  }
}

void EZStereocenter::propagateVertexRemoval(const AtomIndexType removedIndex) {
  auto updateIndex = [&removedIndex](AtomIndexType& index) -> void {
    if(index > removedIndex) {
      --index;
    } else if(index == removedIndex) {
      index = Stereocenter::removalPlaceholder;
    }
  };

  auto newIndex = [&removedIndex](const AtomIndexType& index) -> AtomIndexType {
    assert(index != removedIndex);

    if(index > removedIndex) {
      return index - 1;
    }

    if(index == removedIndex) {
      return Stereocenter::removalPlaceholder;
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

    for(auto& link : ranking.links) {
      link.indexPair = {
        std::min(newIndex(link.indexPair.first), newIndex(link.indexPair.second)),
        std::max(newIndex(link.indexPair.first), newIndex(link.indexPair.second))
      };

      link.cycleSequence = temple::map(
        link.cycleSequence,
        newIndex
      );
    }
  };

  updateRanking(_leftRanking);
  updateRanking(_rightRanking);
}

void EZStereocenter::removeSubstituent(
  const AtomIndexType center,
  const AtomIndexType which
) {
  assert(center == _leftCenter || center == _rightCenter);

  auto dropWhich = [&which](RankingInformation& ranking) {
    for(auto& equalPrioritySet : ranking.sortedSubstituents) {
      temple::inplaceRemoveIf(
        equalPrioritySet,
        [&which](const auto& index) -> bool {
          return index == which;
        }
      );
    }

    // There cannot be any linked pairs for a single substituent
    ranking.links = {};
  };

  if(center == _leftCenter) {
    // Drop which from that ranking
    dropWhich(_leftRanking);
  } else if(center == _rightCenter) {
    dropWhich(_rightRanking);
  }
}


/* Information */
boost::optional<unsigned> EZStereocenter::assigned() const {
  if(_isZOption) {
    return static_cast<unsigned>(
      _isZOption.value()
    );
  }

  return {};
}

boost::optional<unsigned> EZStereocenter::indexOfPermutation() const {
  return assigned();
}

unsigned EZStereocenter::numAssignments() const {
  return numStereopermutations();
}

unsigned EZStereocenter::numStereopermutations() const {
  /* Determine whether there can be two assignments or not. There is only one
   * assignment in the case that on either side, there are two equal
   * substituents
   */
  if(
    (
      _numIndices(_leftRanking) == 2
      && _leftRanking.ligandsRanking.size() == 1
    ) || (
      _numIndices(_rightRanking) == 2
      && _rightRanking.ligandsRanking.size() == 1
    )
  ) {
    return 1;
  }

  return 2;
}

std::vector<DistanceGeometry::ChiralityConstraint> EZStereocenter::chiralityConstraints() const {
  // Three fixed ChiralityConstraints to enforce six-atom coplanarity

  using LigandSequence = DistanceGeometry::ChiralityConstraint::LigandSequence;

  std::vector<DistanceGeometry::ChiralityConstraint> constraints {
    {
      LigandSequence {{
        {_leftHighPriority()},
        {_leftCenter},
        {_rightCenter},
        {_rightHighPriority()}
      }},
      -chiralityConstraintTolerance,
      chiralityConstraintTolerance
    }
  };

  if(_numIndices(_leftRanking) == 2) {
    constraints.emplace_back(
      LigandSequence {{
        {_leftHighPriority()},
        {_leftCenter},
        {_rightCenter},
        {_leftLowPriority()}
      }},
      -chiralityConstraintTolerance,
      chiralityConstraintTolerance
    );
  }

  if(_numIndices(_rightRanking) == 2) {
    constraints.emplace_back(
      LigandSequence {{
        {_rightHighPriority()},
        {_rightCenter},
        {_leftCenter},
        {_rightLowPriority()}
      }},
      -chiralityConstraintTolerance,
      chiralityConstraintTolerance
    );
  }

  return constraints;
}

std::vector<EZStereocenter::DihedralLimits> EZStereocenter::_dihedralLimits(
  const std::function<double(const AtomIndexType)> cycleMultiplierForIndex,
  const double looseningMultiplier
) const {
  /* Whether we have _numStereopermutations == 1 or 2 can be completely ignored because
   * when _numStereopermutations == 1, it is irrelevant which state _isZOption is in,
   * so long as it is set (which it is, already in the constructor). Since both
   * possibilities of considering high-low on one side are equivalent, we just
   * consider the first in all cases. As long as the ranking algorithm works
   * as it should, there should be no issues.
   *
   * So we just make this dependent on the current _isZOption settings.
   */

  std::vector<DihedralLimits> limits;

  auto varianceForSequence = [&](const std::array<AtomIndexType, 4>& sequence) -> double {
    return (
      DistanceGeometry::MoleculeSpatialModel::dihedralAbsoluteVariance
      * looseningMultiplier
      * cycleMultiplierForIndex(sequence.front())
      * cycleMultiplierForIndex(sequence.back())
    );
  };

  // EZStereocenters can impose dihedral limits
  if(_isZOption && _isZOption.value()) {
    // In Z, equal priorities are cis
    for(const auto& dihedralSequence : _equalPriorityDihedralSequences()) {
      limits.emplace_back(
        dihedralSequence,
        std::pair<double, double> {
          0,
          std::min(M_PI, varianceForSequence(dihedralSequence))
        } // cis limit
      );
    }

    for(const auto& dihedralSequence : _differentPriorityDihedralSequences()) {
      limits.emplace_back(
        dihedralSequence,
        std::pair<double, double> {
          std::max(0.0, M_PI - varianceForSequence(dihedralSequence)),
          M_PI
        } // trans
      );
    }
  }

  if(_isZOption && !_isZOption.value()) {
    // In E, equal priorities are trans
    for(const auto& dihedralSequence : _equalPriorityDihedralSequences()) {
      limits.emplace_back(
        dihedralSequence,
        std::pair<double, double> {
          std::max(0.0, M_PI - varianceForSequence(dihedralSequence)),
          M_PI
        } // trans
      );
    }

    for(const auto& dihedralSequence : _differentPriorityDihedralSequences()) {
      limits.emplace_back(
        dihedralSequence,
        std::pair<double, double> {
          0,
          std::min(M_PI, varianceForSequence(dihedralSequence))
        } // cis limit
      );
    }
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

  return limits;
}

const RankingInformation& EZStereocenter::getLeftRanking() const {
  return _leftRanking;
}
const RankingInformation& EZStereocenter::getRightRanking() const {
  return _rightRanking;
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

  if(numStereopermutations() == 1) {
    returnString += ". Is non-stereogenic.";
  } else if(numStereopermutations() == 2) {
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
    "EZ-"s + std::to_string(numStereopermutations())
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

void EZStereocenter::setModelInformation(
  DistanceGeometry::MoleculeSpatialModel& model,
  const std::function<double(const AtomIndexType)> cycleMultiplierForIndex,
  const double looseningMultiplier
) const {
  using ModelType = DistanceGeometry::MoleculeSpatialModel;

  /* Angles */
  auto addAngle = [&](
    const AtomIndexType i,
    const AtomIndexType j,
    const AtomIndexType k
  ) -> void {
    double variance = (
      ModelType::angleAbsoluteVariance
      * looseningMultiplier
      * cycleMultiplierForIndex(i)
      * cycleMultiplierForIndex(k)
    );

    model.setAngleBoundsIfEmpty(
      {{i, j, k}},
      DistanceGeometry::ValueBounds {
        std::max(0.0, 2 * M_PI / 3 - variance),
        std::min(M_PI, 2 * M_PI / 3 + variance)
      }
    );
  };

  // Left
  addAngle(_leftHighPriority(), _leftCenter, _rightCenter);
  if(_numIndices(_leftRanking) == 2) {
    addAngle(_leftLowPriority(), _leftCenter, _rightCenter);
  }

  // Right
  addAngle(_leftCenter, _rightCenter, _rightHighPriority());
  if(_numIndices(_rightRanking) == 2) {
    addAngle(_leftCenter, _rightCenter, _rightLowPriority());
  }

  /* Dihedrals */
  for(
    const auto& dihedralLimit :
    _dihedralLimits(cycleMultiplierForIndex, looseningMultiplier)
  ) {
    model.setDihedralBoundsIfEmpty(
      dihedralLimit.indices,
      dihedralLimit.lower,
      dihedralLimit.upper
    );
  }
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
    && _leftRanking.links == other._leftRanking.links
    && _rightRanking.links == other._rightRanking.links
    && _isZOption == other._isZOption
  );
}

bool EZStereocenter::operator < (const EZStereocenter& other) const {
  return temple::consecutiveCompareSmaller(
    _leftCenter,
    other._leftCenter,
    _rightCenter,
    other._rightCenter,
    _leftRanking.sortedSubstituents,
    other._leftRanking.sortedSubstituents,
    _rightRanking.sortedSubstituents,
    other._rightRanking.sortedSubstituents,
    _isZOption,
    other._isZOption
  );
}

} // namespace Stereocenters

} // namespace molassembler
