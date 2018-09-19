// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/BondStereocenter.h"

#include "stereopermutation/Composites.h"
#include "temple/Adaptors/AllPairs.h"
#include "temple/constexpr/Math.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Functional.h"
#include "temple/OrderedPair.h"
#include "temple/Random.h"

#include "molassembler/AtomStereocenter.h"
#include "molassembler/AngstromWrapper.h"
#include "molassembler/Detail/DelibHelpers.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"
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

constexpr double BondStereocenter::chiralityConstraintTolerance;
constexpr double BondStereocenter::assignmentAcceptanceDihedralThreshold;

/* Impl */
struct BondStereocenter::Impl {
  Impl(
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB,
    BondIndex edge
  );

  void assign(boost::optional<unsigned> assignment);

  void assignRandom();

  void fit(
    const AngstromWrapper& angstromWrapper,
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB
  );

/* Information */
  boost::optional<unsigned> assigned() const;

  bool hasSameCompositeOrientation(const BondStereocenter::Impl& other) const;

  boost::optional<unsigned> indexOfPermutation() const;

  unsigned numAssignments() const;

  unsigned numStereopermutations() const;

  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints(
    double looseningMultiplier,
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB
  ) const;

  std::string info() const;

  std::string rankInfo() const;

  BondIndex edge() const;

  void setModelInformation(
    DistanceGeometry::SpatialModel& model,
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB,
    double looseningMultiplier
  ) const;

/* Operators */
  bool operator == (const Impl& other) const;
  bool operator != (const Impl& other) const;

private:
  stereopermutation::Composite _composite;
  BondIndex _edge;
  boost::optional<unsigned> _assignment;

  //! Yields abstract ligand characters at their symmetry positions
  static std::vector<char> _charifyRankedLigands(
    const std::vector<std::vector<unsigned>>& ligandsRanking,
    const std::vector<unsigned>& symmetryPositionMap
  );

  static stereopermutation::Composite::OrientationState _makeOrientationState(
    const AtomStereocenter& focalStereocenter,
    const AtomStereocenter& attachedStereocenter
  );
};

/* Static Impl member functions */
std::vector<char> BondStereocenter::Impl::_charifyRankedLigands(
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

stereopermutation::Composite::OrientationState BondStereocenter::Impl::_makeOrientationState(
  const AtomStereocenter& focalStereocenter,
  const AtomStereocenter& attachedStereocenter
) {
  return {
    focalStereocenter.getSymmetry(),
    SymmetryMapHelper::getSymmetryPositionOf(
      focalStereocenter.getRanking().getLigandIndexOf(
        attachedStereocenter.centralIndex()
      ),
      focalStereocenter.getSymmetryPositionMap()
    ),
    _charifyRankedLigands(
      focalStereocenter.getRanking().ligandsRanking,
      focalStereocenter.getSymmetryPositionMap()
    ),
    focalStereocenter.centralIndex()
  };
}

/* Constructors */
BondStereocenter::Impl::Impl(
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB,
  const BondIndex edge
) : _composite {
      _makeOrientationState(stereocenterA, stereocenterB),
      _makeOrientationState(stereocenterB, stereocenterA)
    },
    _edge(edge),
    _assignment(boost::none)
{
  assert(
    !stereocenterA.getRanking().hasHapticLigands()
    && !stereocenterB.getRanking().hasHapticLigands()
  );
}

/* Public members */
/* Modification */
void BondStereocenter::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment) {
    assert(assignment.value() < _composite.permutations());
  }

  _assignment = assignment;
}

void BondStereocenter::Impl::assignRandom() {
  assert(_composite.permutations() > 0);

  _assignment = temple::random::getSingle<unsigned>(
    0,
    _composite.permutations() - 1,
    randomnessEngine()
  );
}

void BondStereocenter::Impl::fit(
  const AngstromWrapper& angstromWrapper,
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB
) {
  // Early exit
  if(_composite.permutations() == 0) {
    _assignment = boost::none;
    return;
  }

  // Can the following selections be done with a single branch?
  const AtomStereocenter& firstStereocenter = (
    stereocenterA.centralIndex() == _composite.orientations().first.identifier
    ? stereocenterA
    : stereocenterB
  );

  const AtomStereocenter& secondStereocenter = (
    stereocenterB.centralIndex() == _composite.orientations().second.identifier
    ? stereocenterB
    : stereocenterA
  );

  // For all atoms making up a ligand, decide on the spatial average position
  const std::vector<Eigen::Vector3d> firstLigandPositions = temple::map(
    firstStereocenter.getRanking().ligands,
    [&angstromWrapper](const std::vector<AtomIndex>& ligandAtoms) -> Eigen::Vector3d {
      return DelibHelpers::averagePosition(angstromWrapper.positions, ligandAtoms);
    }
  );

  const std::vector<Eigen::Vector3d> secondLigandPositions = temple::map(
    secondStereocenter.getRanking().ligands,
    [&angstromWrapper](const std::vector<AtomIndex>& ligandAtoms) -> Eigen::Vector3d {
      return DelibHelpers::averagePosition(angstromWrapper.positions, ligandAtoms);
    }
  );

  unsigned firstSymmetryPosition, secondSymmetryPosition;
  double dihedralAngle;

  double bestPenalty = std::numeric_limits<double>::max();
  std::vector<unsigned> bestAssignment;

  /* TODO long-term, as soon as this is proven to work as it should, this
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
        firstStereocenter.getSymmetryPositionMap()
      );
      unsigned secondLigandIndex = SymmetryMapHelper::getLigandIndexAt(
        secondSymmetryPosition,
        secondStereocenter.getSymmetryPositionMap()
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
        angstromWrapper.positions.at(firstStereocenter.centralIndex()).toEigenVector(),
        angstromWrapper.positions.at(secondStereocenter.centralIndex()).toEigenVector(),
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
boost::optional<unsigned> BondStereocenter::Impl::assigned() const {
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

bool BondStereocenter::Impl::hasSameCompositeOrientation(const BondStereocenter::Impl& other) const {
  return _composite == other._composite;
}

boost::optional<unsigned> BondStereocenter::Impl::indexOfPermutation() const {
  if(_assignment && _composite.isIsotropic()) {
    return 0u;
  }

  return _assignment;
}

unsigned BondStereocenter::Impl::numAssignments() const {
  if(_composite.isIsotropic()) {
    return 1;
  }

  return _composite.permutations();
}

unsigned BondStereocenter::Impl::numStereopermutations() const {
  if(_composite.isIsotropic()) {
    return 1;
  }

  return _composite.permutations();
}

std::vector<DistanceGeometry::ChiralityConstraint> BondStereocenter::Impl::chiralityConstraints(
  const double looseningMultiplier,
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB
) const {
  // Can the following selections be done with a single branch?
  const AtomStereocenter& firstStereocenter = (
    stereocenterA.centralIndex() == _composite.orientations().first.identifier
    ? stereocenterA
    : stereocenterB
  );

  const AtomStereocenter& secondStereocenter = (
    stereocenterB.centralIndex() == _composite.orientations().second.identifier
    ? stereocenterB
    : stereocenterA
  );

  std::vector<DistanceGeometry::ChiralityConstraint> constraints;

  unsigned firstSymmetryPosition, secondSymmetryPosition;
  double dihedralAngle;

  for(const auto& dihedralTuple : _composite.dihedrals(_assignment.value())) {
    std::tie(firstSymmetryPosition, secondSymmetryPosition, dihedralAngle) = dihedralTuple;

    if(std::fabs(dihedralAngle) < M_PI / 180.0 || std::fabs(M_PI - std::fabs(dihedralAngle)) < M_PI / 180.0) {
      constraints.emplace_back(
        DistanceGeometry::ChiralityConstraint::LigandSequence {
          firstStereocenter.getRanking().ligands.at(
            SymmetryMapHelper::getLigandIndexAt(
              firstSymmetryPosition,
              firstStereocenter.getSymmetryPositionMap()
            )
          ),
          {firstStereocenter.centralIndex()},
          {secondStereocenter.centralIndex()},
          secondStereocenter.getRanking().ligands.at(
            SymmetryMapHelper::getLigandIndexAt(
              secondSymmetryPosition,
              secondStereocenter.getSymmetryPositionMap()
            )
          )
        },
        - chiralityConstraintTolerance * looseningMultiplier,
        chiralityConstraintTolerance * looseningMultiplier
      );
    }
  }

  return constraints;
}

std::string BondStereocenter::Impl::info() const {
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

std::string BondStereocenter::Impl::rankInfo() const {
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

BondIndex BondStereocenter::Impl::edge() const {
  // Return a standard form of smaller first
  return _edge;
}

void BondStereocenter::Impl::setModelInformation(
  DistanceGeometry::SpatialModel& model,
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB,
  const double looseningMultiplier
) const {
  // Can the following selections be done with a single branch?
  const AtomStereocenter& firstStereocenter = (
    stereocenterA.centralIndex() == _composite.orientations().first.identifier
    ? stereocenterA
    : stereocenterB
  );

  const AtomStereocenter& secondStereocenter = (
    stereocenterB.centralIndex() == _composite.orientations().second.identifier
    ? stereocenterB
    : stereocenterA
  );

  using ModelType = DistanceGeometry::SpatialModel;

  assert(_assignment);

  /* Only dihedrals */
  unsigned firstSymmetryPosition, secondSymmetryPosition;
  double dihedralAngle;

  for(const auto& dihedralTuple : _composite.dihedrals(_assignment.value())) {
    std::tie(firstSymmetryPosition, secondSymmetryPosition, dihedralAngle) = dihedralTuple;

    unsigned firstLigandIndex = SymmetryMapHelper::getLigandIndexAt(
      firstSymmetryPosition,
      firstStereocenter.getSymmetryPositionMap()
    );
    unsigned secondLigandIndex = SymmetryMapHelper::getLigandIndexAt(
      secondSymmetryPosition,
      secondStereocenter.getSymmetryPositionMap()
    );

    /* A simple central plus-minus variance calculation for the dihedral would
     * yield values outside the SpatialModel dihedralClampBounds, which are [0,
     * M_PI] (since d(phi) is y-symmetric), so we have some adjusting to do.
     */
    auto reduceToBounds = [](double phi) -> double {
      // Perform 2pi adjustments until we are in (-pi, pi]
      if(phi <= -M_PI) {
        unsigned n = std::floor((M_PI + std::fabs(phi)) / (2 * M_PI));
        phi += n * 2 * M_PI;
      }

      if(phi > M_PI) {
        unsigned n = std::floor((M_PI + phi) / (2 * M_PI));
        phi -= n * 2 * M_PI;
      }

      // Within (-pi, pi], d(phi) is y-axis-symmetric
      return std::fabs(phi);
    };

    std::array<double, 3> possibleBoundaryValues {
      {
        reduceToBounds(dihedralAngle - ModelType::dihedralAbsoluteVariance * looseningMultiplier),
        reduceToBounds(dihedralAngle),
        reduceToBounds(dihedralAngle + ModelType::dihedralAbsoluteVariance * looseningMultiplier)
      }
    };

    // The constructor will swap passed values so that the smaller becomes .first
    temple::OrderedPair<double> boundedDihedralValues {
      *std::min_element(std::begin(possibleBoundaryValues), std::end(possibleBoundaryValues)),
      *std::max_element(std::begin(possibleBoundaryValues), std::end(possibleBoundaryValues))
    };

    // Take the ordered values
    DistanceGeometry::ValueBounds dihedralBounds {
      boundedDihedralValues.first,
      boundedDihedralValues.second
    };

    temple::forEach(
      temple::adaptors::allPairs(
        firstStereocenter.getRanking().ligands.at(firstLigandIndex),
        secondStereocenter.getRanking().ligands.at(secondLigandIndex)
      ),
      [&](const AtomIndex firstIndex, const AtomIndex secondIndex) -> void {
        model.setDihedralBoundsIfEmpty(
          std::array<AtomIndex, 4> {
            firstIndex,
            firstStereocenter.centralIndex(),
            secondStereocenter.centralIndex(),
            secondIndex
          },
          dihedralBounds
        );
      }
    );
  }
}

/* Operators */
bool BondStereocenter::Impl::operator == (const Impl& other) const {
  return (
    _composite == other._composite
    && _assignment == other._assignment
  );
}

bool BondStereocenter::Impl::operator != (const Impl& other) const {
  return !(*this == other);
}


/* BondStereocenter implementation */
BondStereocenter::BondStereocenter(BondStereocenter&& other) noexcept = default;
BondStereocenter& BondStereocenter::operator = (BondStereocenter&& other) noexcept = default;

BondStereocenter::BondStereocenter(const BondStereocenter& other) : _pImpl(
  std::make_unique<Impl>(*other._pImpl)
) {}
BondStereocenter& BondStereocenter::operator = (const BondStereocenter& other) {
  *_pImpl = *other._pImpl;
  return *this;
}
BondStereocenter::~BondStereocenter() = default;

BondStereocenter::BondStereocenter(
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB,
  const BondIndex& edge
) {
  _pImpl = std::make_unique<Impl>(
    stereocenterA,
    stereocenterB,
    edge
  );
}

void BondStereocenter::assign(boost::optional<unsigned> assignment) {
  _pImpl -> assign(std::move(assignment));
}

void BondStereocenter::assignRandom() {
  _pImpl -> assignRandom();
}

void BondStereocenter::fit(
  const AngstromWrapper& angstromWrapper,
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB
) {
  _pImpl -> fit(
    angstromWrapper,
    stereocenterA,
    stereocenterB
  );
}

boost::optional<unsigned> BondStereocenter::assigned() const {
  return _pImpl -> assigned();
}

bool BondStereocenter::hasSameCompositeOrientation(const BondStereocenter& other) const {
  return _pImpl -> hasSameCompositeOrientation(*other._pImpl);
}

boost::optional<unsigned> BondStereocenter::indexOfPermutation() const {
  return _pImpl -> indexOfPermutation();
}

unsigned BondStereocenter::numAssignments() const {
  return _pImpl -> numAssignments();
}

unsigned BondStereocenter::numStereopermutations() const {
  return _pImpl -> numStereopermutations();
}

std::vector<DistanceGeometry::ChiralityConstraint> BondStereocenter::chiralityConstraints(
  double looseningMultiplier,
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB
) const {
  return _pImpl -> chiralityConstraints(
    looseningMultiplier,
    stereocenterA,
    stereocenterB
  );
}

std::string BondStereocenter::info() const {
  return _pImpl -> info();
}

std::string BondStereocenter::rankInfo() const {
  return _pImpl -> rankInfo();
}

BondIndex BondStereocenter::edge() const {
  return _pImpl -> edge();
}

void BondStereocenter::setModelInformation(
  DistanceGeometry::SpatialModel& model,
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB,
  double looseningMultiplier
) const {
  _pImpl -> setModelInformation(model, stereocenterA, stereocenterB, looseningMultiplier);
}

bool BondStereocenter::operator == (const BondStereocenter& other) const {
  return *_pImpl == *other._pImpl;
}

bool BondStereocenter::operator != (const BondStereocenter& other) const {
  return !(
    *_pImpl == *other._pImpl
  );
}

} // namespace molassembler
