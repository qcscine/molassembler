#include "molassembler/Stereopermutators/BondStereopermutatorImpl.h"

#include "temple/Adaptors/AllPairs.h"
#include "temple/constexpr/Math.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Functional.h"
#include "temple/OrderedPair.h"
#include "temple/Random.h"

#include "molassembler/AtomStereopermutator.h"
#include "molassembler/AngstromWrapper.h"
#include "molassembler/Detail/DelibHelpers.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/Molecule.h"
#include "molassembler/RankingInformation.h"

namespace molassembler {

struct SymmetryMapHelper {
  static unsigned getSymmetryPositionOf(unsigned ligandIndex, const std::vector<unsigned>& map) {
    return map.at(ligandIndex);
  }

  static unsigned getLigandIndexAt(unsigned symmetryPosition, const std::vector<unsigned>& map) {
    auto findIter = std::find(
      std::begin(map),
      std::end(map),
      symmetryPosition
    );

    assert(findIter != std::end(map));

    return findIter - std::begin(map);
  }
};

std::vector<char> BondStereopermutator::Impl::_charifyRankedLigands(
  const std::vector<std::vector<unsigned>>& ligandsRanking,
  const std::vector<unsigned>& symmetryPositionMap
) {
  std::vector<char> characters (symmetryPositionMap.size());

  char currentChar = 'A';
  for(const auto& equalPrioritySet : ligandsRanking) {
    for(const auto& ligandIndex : equalPrioritySet) {
      characters.at(
        symmetryPositionMap.at(ligandIndex)
      ) = currentChar;
    }

    ++currentChar;
  }

  return characters;
}

const stereopermutation::Composite& BondStereopermutator::Impl::composite() const {
  return _composite;
}

stereopermutation::Composite::OrientationState BondStereopermutator::Impl::_makeOrientationState(
  const AtomStereopermutator& focalStereopermutator,
  const AtomStereopermutator& attachedStereopermutator
) {
  return {
    focalStereopermutator.getSymmetry(),
    SymmetryMapHelper::getSymmetryPositionOf(
      focalStereopermutator.getRanking().getLigandIndexOf(
        attachedStereopermutator.centralIndex()
      ),
      focalStereopermutator.getSymmetryPositionMap()
    ),
    _charifyRankedLigands(
      focalStereopermutator.getRanking().ligandsRanking,
      focalStereopermutator.getSymmetryPositionMap()
    ),
    focalStereopermutator.centralIndex()
  };
}

/* Constructors */
BondStereopermutator::Impl::Impl(
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const BondIndex edge
) : _composite {
      _makeOrientationState(stereopermutatorA, stereopermutatorB),
      _makeOrientationState(stereopermutatorB, stereopermutatorA)
    },
    _edge(edge),
    _assignment(boost::none)
{
  assert(
    !stereopermutatorA.getRanking().hasHapticLigands()
    && !stereopermutatorB.getRanking().hasHapticLigands()
  );
}

/* Public members */
/* Modification */
void BondStereopermutator::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment) {
    assert(assignment.value() < _composite.permutations());
  }

  _assignment = assignment;
}

void BondStereopermutator::Impl::assignRandom() {
  assert(_composite.permutations() > 0);

  _assignment = temple::random::getSingle<unsigned>(
    0,
    _composite.permutations() - 1,
    randomnessEngine()
  );
}

void BondStereopermutator::Impl::fit(
  const AngstromWrapper& angstromWrapper,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB
) {
  // Early exit
  if(_composite.permutations() == 0) {
    _assignment = boost::none;
    return;
  }

  // Can the following selections be done with a single branch?
  const AtomStereopermutator& firstStereopermutator = (
    stereopermutatorA.centralIndex() == _composite.orientations().first.identifier
    ? stereopermutatorA
    : stereopermutatorB
  );

  const AtomStereopermutator& secondStereopermutator = (
    stereopermutatorB.centralIndex() == _composite.orientations().second.identifier
    ? stereopermutatorB
    : stereopermutatorA
  );

  // For all atoms making up a ligand, decide on the spatial average position
  const std::vector<Eigen::Vector3d> firstLigandPositions = temple::map(
    firstStereopermutator.getRanking().ligands,
    [&angstromWrapper](const std::vector<AtomIndex>& ligandAtoms) -> Eigen::Vector3d {
      return DelibHelpers::averagePosition(angstromWrapper.positions, ligandAtoms);
    }
  );

  const std::vector<Eigen::Vector3d> secondLigandPositions = temple::map(
    secondStereopermutator.getRanking().ligands,
    [&angstromWrapper](const std::vector<AtomIndex>& ligandAtoms) -> Eigen::Vector3d {
      return DelibHelpers::averagePosition(angstromWrapper.positions, ligandAtoms);
    }
  );

  unsigned firstSymmetryPosition, secondSymmetryPosition;
  double dihedralAngle;

  double bestPenalty = std::numeric_limits<double>::max();
  std::vector<unsigned> bestAssignment;

  /*! @todo long-term, as soon as this is proven to work as it should, this
   * computation can be abbreviated to exit the dihedralTuple loop as soon as
   * the threshold is surpassed to avoid unnecessary dihedral computations and
   * domain corrections.
   */

  for(unsigned i = 0; i < _composite.permutations(); ++i) {
    double penalty = 0.0;
    for(const auto& dihedralTuple : _composite.dihedrals(i)) {
      std::tie(firstSymmetryPosition, secondSymmetryPosition, dihedralAngle) = dihedralTuple;

      // Get ligand index of leftSymmetryPosition in left
      unsigned firstLigandIndex = SymmetryMapHelper::getLigandIndexAt(
        firstSymmetryPosition,
        firstStereopermutator.getSymmetryPositionMap()
      );
      unsigned secondLigandIndex = SymmetryMapHelper::getLigandIndexAt(
        secondSymmetryPosition,
        secondStereopermutator.getSymmetryPositionMap()
      );

      /* Dihedral angle differences aren't as easy as |b - a|, since
       * dihedrals are defined over (-pi, pi], so in the worst case
       *
       *   a = pi, b = -pi + delta
       *   |(-pi + delta) - pi| = 2 pi - delta.
       *
       * Any difference has to be modified by a 2 pi-periodic function
       * calculating the distance to the nearest multiple of pi (like a
       * triangle wave function) or we have to ensure that the difference is
       * also in the definition interval and then calculate the distance to
       * zero.
       */

      double dihedralDifference = DelibHelpers::dihedral(
        firstLigandPositions.at(firstLigandIndex),
        angstromWrapper.positions.at(firstStereopermutator.centralIndex()).toEigenVector(),
        angstromWrapper.positions.at(secondStereopermutator.centralIndex()).toEigenVector(),
        secondLigandPositions.at(secondLigandIndex)
      ) - dihedralAngle;

      // + pi is part of the definition interval, so use greater than
      if(dihedralDifference > M_PI) {
        dihedralDifference -= 2 * M_PI;
      }

      // - pi is not part of the definition interval, so use less than or equal
      if(dihedralDifference <= -M_PI) {
        dihedralDifference += 2 * M_PI;
      }

      penalty += std::fabs(dihedralDifference);
    }

    // Check if this penalty is within acceptance threshold
    if(penalty > assignmentAcceptanceDihedralThreshold * _composite.dihedrals(i).size()) {
      continue;
    }

    if(penalty < bestPenalty) {
      bestPenalty = penalty;
      bestAssignment = {i};
    } else if(penalty == bestPenalty) {
      bestAssignment.push_back(i);
    }
  }

  if(bestAssignment.size() == 1) {
    _assignment = bestAssignment.front();
  }
}

/* Information */
boost::optional<unsigned> BondStereopermutator::Impl::assigned() const {
  /* If the underlying composite is isotropic, it does not matter which of those
   * permutations by symmetry position is the factual spatial arrangement (since
   * they are all rotationally equivalent). We have to spoof that there is only
   * one arrangement in this case (although we need all of them for spatial
   * fitting).
   */
  if(_assignment && _composite.isIsotropic()) {
    return 0u;
  }

  return _assignment;
}

bool BondStereopermutator::Impl::hasSameCompositeOrientation(const BondStereopermutator::Impl& other) const {
  return _composite == other._composite;
}

boost::optional<unsigned> BondStereopermutator::Impl::indexOfPermutation() const {
  if(_assignment && _composite.isIsotropic()) {
    return 0u;
  }

  return _assignment;
}

unsigned BondStereopermutator::Impl::numAssignments() const {
  if(_composite.isIsotropic()) {
    return 1;
  }

  return _composite.permutations();
}

unsigned BondStereopermutator::Impl::numStereopermutations() const {
  if(_composite.isIsotropic()) {
    return 1;
  }

  return _composite.permutations();
}

std::vector<DistanceGeometry::ChiralityConstraint> BondStereopermutator::Impl::chiralityConstraints(
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB
) const {
  // Can the following selections be done with a single branch?
  const AtomStereopermutator& firstStereopermutator = (
    stereopermutatorA.centralIndex() == _composite.orientations().first.identifier
    ? stereopermutatorA
    : stereopermutatorB
  );

  const AtomStereopermutator& secondStereopermutator = (
    stereopermutatorB.centralIndex() == _composite.orientations().second.identifier
    ? stereopermutatorB
    : stereopermutatorA
  );

  std::vector<DistanceGeometry::ChiralityConstraint> constraints;

  unsigned firstSymmetryPosition, secondSymmetryPosition;
  double dihedralAngle;

  for(const auto& dihedralTuple : _composite.dihedrals(_assignment.value())) {
    std::tie(firstSymmetryPosition, secondSymmetryPosition, dihedralAngle) = dihedralTuple;

    if(std::fabs(dihedralAngle) < M_PI / 180.0 || std::fabs(M_PI - std::fabs(dihedralAngle)) < M_PI / 180.0) {
      constraints.emplace_back(
        DistanceGeometry::ChiralityConstraint::LigandSequence {
          firstStereopermutator.getRanking().ligands.at(
            SymmetryMapHelper::getLigandIndexAt(
              firstSymmetryPosition,
              firstStereopermutator.getSymmetryPositionMap()
            )
          ),
          {firstStereopermutator.centralIndex()},
          {secondStereopermutator.centralIndex()},
          secondStereopermutator.getRanking().ligands.at(
            SymmetryMapHelper::getLigandIndexAt(
              secondSymmetryPosition,
              secondStereopermutator.getSymmetryPositionMap()
            )
          )
        },
        0,
        0
      );
    }
  }

  return constraints;
}

std::string BondStereopermutator::Impl::info() const {
  using namespace std::string_literals;

  std::string returnString =  "B on "s;

  returnString += std::to_string(_composite.orientations().first.identifier);
  returnString += "-";
  returnString += std::to_string(_composite.orientations().second.identifier);

  if(numStereopermutations() == 1) {
    returnString += ". Is non-stereogenic.";
  } else {
    returnString += ". Is "s;
    if(_assignment) {
      returnString += std::to_string(_assignment.value());
    } else {
      returnString += "u";
    }

    returnString += "/"s + std::to_string(numStereopermutations());
  }

  return returnString;
}

std::string BondStereopermutator::Impl::rankInfo() const {
  using namespace std::string_literals;

  return (
    "B-"s + std::to_string(numStereopermutations())
    + "-"s + (
      assigned()
      ? std::to_string(assigned().value())
      : "u"s
    )
  );
}

BondIndex BondStereopermutator::Impl::edge() const {
  // Return a standard form of smaller first
  return _edge;
}

/* Operators */
bool BondStereopermutator::Impl::operator == (const Impl& other) const {
  return (
    _composite == other._composite
    && _assignment == other._assignment
  );
}

bool BondStereopermutator::Impl::operator != (const Impl& other) const {
  return !(*this == other);
}


} // namespace molassembler
