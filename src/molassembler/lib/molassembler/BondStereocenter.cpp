#include "temple/constexpr/ConsecutiveCompare.h"
#include "temple/constexpr/Math.h"
#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Random.h"

#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/SpatialModel.h"

#include "detail/DelibHelpers.h"
#include "AtomStereocenter.h"
#include "BondStereocenter.h"
#include "BondDistance.h"
#include "Molecule.h"
#include "RankingInformation.h"

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

bool BondStereocenter::AtomState::operator == (const BondStereocenter::AtomState& other) const {
  return (
    index == other.index
    && symmetry == other.symmetry
  );
}

bool BondStereocenter::AtomState::operator != (const BondStereocenter::AtomState& other) const {
  return !(*this == other);
}

std::vector<char> BondStereocenter::_charifyRankedLigands(const std::vector<std::vector<unsigned>> ligandsRanking) {
  std::vector<char> characters;

  char currentChar = 'A';
  for(const auto& equalPrioritySet : ligandsRanking) {
    for(unsigned i = 0; i < equalPrioritySet.size(); ++i) {
      characters.push_back(currentChar);
    }

    ++currentChar;
  }

  return characters;
}

/* Constructors */
BondStereocenter::BondStereocenter(
  const AtomStereocenter& leftStereocenter,
  const AtomStereocenter& rightStereocenter,
  const GraphType::edge_descriptor edge
) : _composite {
      leftStereocenter.getSymmetry(),
      rightStereocenter.getSymmetry(),
      SymmetryMapHelper::getSymmetryPositionOf(
        leftStereocenter.getRanking().getLigandIndexOf(
          rightStereocenter.centralIndex()
        ),
        leftStereocenter.getSymmetryPositionMap()
      ),
      SymmetryMapHelper::getSymmetryPositionOf(
        rightStereocenter.getRanking().getLigandIndexOf(
          leftStereocenter.centralIndex()
        ),
        rightStereocenter.getSymmetryPositionMap()
      ),
      _charifyRankedLigands(leftStereocenter.getRanking().ligandsRanking),
      _charifyRankedLigands(rightStereocenter.getRanking().ligandsRanking)
    },
    _edge(edge),
    _assignment(boost::none)
{
  assert(
    !leftStereocenter.getRanking().hasHapticLigands()
    && !rightStereocenter.getRanking().hasHapticLigands()
  );

  // Copy AtomStereocenter state so we can handle state propagation eventually
  _left.index = leftStereocenter.centralIndex();
  _left.symmetry = leftStereocenter.getSymmetry();
  // _left.ranking = leftStereocenter.getRanking();
  _right.index = rightStereocenter.centralIndex();
  _right.symmetry = rightStereocenter.getSymmetry();
  // _right.ranking = rightStereocenter.getRanking();
}

/* Public members */
/* Modification */
void BondStereocenter::assign(boost::optional<unsigned> assignment) {
  if(assignment) {
    assert(assignment.value() < _composite.permutations());
  }

  _assignment = assignment;
}

void BondStereocenter::assignRandom() {
  assert(_composite.permutations() > 0);

  _assignment = prng.getSingle<unsigned>(0, _composite.permutations() - 1);
}

void BondStereocenter::fit(
  const AngstromWrapper& angstromWrapper,
  const AtomStereocenter& left,
  const AtomStereocenter& right
) {
  // For all atoms making up a ligand, decide on the spatial average position
  const std::vector<Eigen::Vector3d> leftLigandPositions = temple::mapToVector(
    left.getRanking().ligands,
    [&angstromWrapper](const std::vector<AtomIndexType>& ligandAtoms) -> Eigen::Vector3d {
      return DelibHelpers::averagePosition(angstromWrapper.positions, ligandAtoms);
    }
  );

  const std::vector<Eigen::Vector3d> rightLigandPositions = temple::mapToVector(
    right.getRanking().ligands,
    [&angstromWrapper](const std::vector<AtomIndexType>& ligandAtoms) -> Eigen::Vector3d {
      return DelibHelpers::averagePosition(angstromWrapper.positions, ligandAtoms);
    }
  );

  unsigned leftSymmetryPosition, rightSymmetryPosition;
  double dihedralAngle;

  std::vector<double> penalties (_composite.permutations(), 0.0);

  for(unsigned i = 0; i < _composite.permutations(); ++i) {
    for(const auto& dihedralTuple : _composite.dihedrals(i)) {
      std::tie(leftSymmetryPosition, rightSymmetryPosition, dihedralAngle) = dihedralTuple;

      // Get ligand index of leftSymmetryPosition in left
      unsigned leftLigandIndex = SymmetryMapHelper::getLigandIndexAt(
        leftSymmetryPosition,
        left.getSymmetryPositionMap()
      );
      unsigned rightLigandIndex = SymmetryMapHelper::getLigandIndexAt(
        rightSymmetryPosition,
        right.getSymmetryPositionMap()
      );

      penalties.at(i) += std::fabs(
        DelibHelpers::dihedral(
          leftLigandPositions.at(leftLigandIndex),
          angstromWrapper.positions.at(left.centralIndex()).toEigenVector(),
          angstromWrapper.positions.at(right.centralIndex()).toEigenVector(),
          rightLigandPositions.at(rightLigandIndex)
        ) - dihedralAngle
      );
    }
  }

  _assignment = std::min_element(
    std::begin(penalties),
    std::end(penalties)
  ) - std::begin(penalties);
}

/* Information */
boost::optional<unsigned> BondStereocenter::assigned() const {
  return _assignment;
}

boost::optional<unsigned> BondStereocenter::indexOfPermutation() const {
  return _assignment;
}

unsigned BondStereocenter::numAssignments() const {
  return _composite.permutations();
}

unsigned BondStereocenter::numStereopermutations() const {
  return _composite.permutations();
}

std::vector<DistanceGeometry::ChiralityConstraint> BondStereocenter::chiralityConstraints(
  const double looseningMultiplier,
  const AtomStereocenter& left,
  const AtomStereocenter& right
) const {
  std::vector<DistanceGeometry::ChiralityConstraint> constraints;

  unsigned leftSymmetryPosition, rightSymmetryPosition;
  double dihedralAngle;

  for(const auto& dihedralTuple : _composite.dihedrals(_assignment.value())) {
    std::tie(leftSymmetryPosition, rightSymmetryPosition, dihedralAngle) = dihedralTuple;

    if(dihedralAngle == 0.0 || dihedralAngle == M_PI) {
      constraints.emplace_back(
        DistanceGeometry::ChiralityConstraint::LigandSequence {
          left.getRanking().ligands.at(
            SymmetryMapHelper::getLigandIndexAt(
              leftSymmetryPosition,
              left.getSymmetryPositionMap()
            )
          ),
          {left.centralIndex()},
          {right.centralIndex()},
          right.getRanking().ligands.at(
            SymmetryMapHelper::getLigandIndexAt(
              rightSymmetryPosition,
              right.getSymmetryPositionMap()
            )
          )
        },
        -chiralityConstraintTolerance * looseningMultiplier,
        chiralityConstraintTolerance * looseningMultiplier
      );
    }
  }

  return constraints;
}

std::string BondStereocenter::info() const {
  using namespace std::string_literals;

  std::string returnString =  "B on "s;

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

std::string BondStereocenter::rankInfo() const {
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

GraphType::edge_descriptor BondStereocenter::edge() const {
  // Return a standard form of smaller first
  return _edge;
}

void BondStereocenter::setModelInformation(
  DistanceGeometry::SpatialModel& model,
  const AtomStereocenter& left,
  const AtomStereocenter& right,
  const double looseningMultiplier
) const {
  using ModelType = DistanceGeometry::SpatialModel;

  assert(_assignment);

  /* Only dihedrals */
  unsigned leftSymmetryPosition, rightSymmetryPosition;
  double dihedralAngle;

  for(const auto& dihedralTuple : _composite.dihedrals(_assignment.value())) {
    std::tie(leftSymmetryPosition, rightSymmetryPosition, dihedralAngle) = dihedralTuple;

    unsigned leftLigandIndex = SymmetryMapHelper::getLigandIndexAt(
      leftSymmetryPosition,
      left.getSymmetryPositionMap()
    );
    unsigned rightLigandIndex = SymmetryMapHelper::getLigandIndexAt(
      rightSymmetryPosition,
      right.getSymmetryPositionMap()
    );

    temple::forAllPairs(
      left.getRanking().ligands.at(leftLigandIndex),
      right.getRanking().ligands.at(rightLigandIndex),
      [&](const AtomIndexType leftIndex, const AtomIndexType rightIndex) -> void {
        model.setDihedralBoundsIfEmpty(
          std::array<AtomIndexType, 4> {
            leftIndex,
            left.centralIndex(),
            right.centralIndex(),
            rightIndex
          },
          ModelType::makeBoundsFromCentralValue(
            dihedralAngle,
            ModelType::dihedralAbsoluteVariance
            * looseningMultiplier
          )
        );
      }
    );
  }
}

/* Operators */
bool BondStereocenter::operator == (const BondStereocenter& other) const {
  return (
    _composite == other._composite
    && _assignment == other._assignment
  );
}

bool BondStereocenter::operator != (const BondStereocenter& other) const {
  return !(*this == other);
}

} // namespace molassembler
