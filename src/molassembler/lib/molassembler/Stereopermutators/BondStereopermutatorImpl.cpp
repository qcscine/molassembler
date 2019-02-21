/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/BondStereopermutatorImpl.h"

#include "chemical_symmetries/Symmetries.h"

#include "temple/Adaptors/AllPairs.h"
#include "temple/constexpr/Math.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Functional.h"
#include "temple/OrderedPair.h"
#include "temple/Random.h"

#include "molassembler/AngstromWrapper.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/Detail/DelibHelpers.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/Molecule.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/Stereopermutators/PermutationState.h"

namespace Scine {

namespace molassembler {

namespace detail {

bool piPeriodicFPCompare(const double a, const double b) {
  /* Reduces everything to [-pi, pi) bounds and then compares
   */
  auto reduceToBounds = [](const double x) -> double {
    return x - std::floor((x + M_PI)/(2 * M_PI)) * 2 * M_PI;
  };

  return std::fabs(reduceToBounds(a) - reduceToBounds(b)) < 1e-10;
}

/*
 * @note If we use tuple_element_t for SFINAE here, this works only for pair and
 * tuple. This way, it works for all types that have a first and second member,
 * of which the relevant types are std::pair and temple::OrderedPair here.
 */
template<typename HomogeneousPair>
auto& select(
  HomogeneousPair& pair,
  bool selectFirst,
  std::enable_if_t<
    std::is_same<
      decltype(std::declval<HomogeneousPair>().first),
      decltype(std::declval<HomogeneousPair>().second)
    >::value,
    int
  >* /* enableIfPtr */ = nullptr
) {
  if(selectFirst) {
    return pair.first;
  }

  return pair.second;
}

template<typename T1, typename T2, typename ... TupleArgs>
const auto& select(
  const std::tuple<T1, T2, TupleArgs...>& tuple,
  bool selectFirst,
  std::enable_if_t<
    std::is_same<T1, T2>::value,
    int
  >* /* enableIfPtr */ = nullptr
) {
  if(selectFirst) {
    return std::get<0>(tuple);
  }

  return std::get<1>(tuple);
}

} // namespace detail

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

stereopermutation::Composite::OrientationState
BondStereopermutator::Impl::_makeOrientationState(
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
  if(
    stereopermutatorA.getRanking().hasHapticLigands()
    || stereopermutatorB.getRanking().hasHapticLigands()
  ) {
    throw std::logic_error("BondStereopermutators do not support haptic ligands yet.");
  }
}

/* Public members */
/* Modification */
void BondStereopermutator::Impl::assign(boost::optional<unsigned> assignment) {
  if(assignment) {
    /* The distinction here between numAssignments and _composite.permutations()
     * is important because if the composite is isotropic, we simulate that
     * there is only a singular assignment (see numAssignments), although
     * all permutations are generated and present for fitting anyway.
     *
     * If this were _composite.permutations(), we would accept assignments other
     * than zero if the underlying composite is isotropic, but not yield the
     * same assignment index when asked which assignment is currently set
     * in assigned().
     */
    assert(assignment.value() < numAssignments());
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

void BondStereopermutator::Impl::applyPermutation(const std::vector<AtomIndex>& permutation) {
  /* Composite's OrientationState identifiers change (these have no impact on
   * assignment ordering however, so that should be safe)
   */
  _composite.applyIdentifierPermutation(permutation);

  // Edge we sit on changes according to the permutation
  _edge = BondIndex {
    permutation.at(_edge.first),
    permutation.at(_edge.second)
  };

  // Assignment does not change
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
        angstromWrapper.positions.row(firstStereopermutator.centralIndex()),
        angstromWrapper.positions.row(secondStereopermutator.centralIndex()),
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

void BondStereopermutator::Impl::propagateGraphChange(
  const AtomStereopermutatorPropagatedState& oldPermutatorState,
  const AtomStereopermutator& newPermutator
) {
  const RankingInformation& oldRanking = std::get<0>(oldPermutatorState);
  const PermutationState& oldPermutationState = std::get<1>(oldPermutatorState);

  // We assume that the supplied permutators (or their state) were assigned
  assert(std::get<2>(oldPermutatorState));
  assert(newPermutator.assigned());

  /* We assume the old and new symmetry are of the same size (i.e. this is
   * a ranking change propagation, not a substituent addition / removal
   * propagation)
   */
  assert(oldRanking.ligands.size() == newPermutator.getRanking().ligands.size());

  using OrientationState = stereopermutation::Composite::OrientationState;

  /* Construct a new Composite with the new information */
  bool changedIsFirstInOldOrientations = (
    _composite.orientations().first.identifier == newPermutator.centralIndex()
  );

  const OrientationState& oldOrientation = detail::select(
    _composite.orientations(),
    changedIsFirstInOldOrientations
  );

  // Reuse the OrientationState of the "other" atom stereopermutator
  const OrientationState& unchangedOrientation = detail::select(
    _composite.orientations(),
    !changedIsFirstInOldOrientations
  );

  // Generate a new OrientationState for the modified stereopermutator
  stereopermutation::Composite::OrientationState possiblyModifiedOrientation {
    newPermutator.getSymmetry(),
    SymmetryMapHelper::getSymmetryPositionOf(
      newPermutator.getRanking().getLigandIndexOf(
        unchangedOrientation.identifier
      ),
      newPermutator.getSymmetryPositionMap()
    ),
    _charifyRankedLigands(
      newPermutator.getRanking().ligandsRanking,
      newPermutator.getSymmetryPositionMap()
    ),
    newPermutator.centralIndex()
  };

  // In case nothing has changed, we are done and can stop
  if(oldOrientation == possiblyModifiedOrientation) {
    return;
  }

  /* Generate a new set of permutations (note this may reorder its
   * arguments into the orientations() member)
   */
  stereopermutation::Composite newComposite {
    possiblyModifiedOrientation,
    unchangedOrientation
  };

  /* If the new composite is isotropic, there is no reason to try to find
   * an assignment, we can choose any.
   */
  if(newComposite.isIsotropic()) {
    _composite = std::move(newComposite);
    _assignment = 0;
    return;
  }

  /* If this BondStereopermutator is unassigned, and the permutator is not
   * isotropic after the ranking change, then this permutator stays unassigned.
   */
  if(_assignment == boost::none) {
    _composite = std::move(newComposite);
    return;
  }

  /* Find the old permutation in the set of new permutations
   * - Since composite only offers dihedral information in terms of symmetry
   *   positions, we have to translate these back into ligand indices
   * - Transform symmetry positions through the sequence
   *
   *   old symmetry position
   *   -> ligand index (using old permutation state)
   *   -> atom index (using old ranking)
   *   -> ligand index (using new ranking)
   *   -> new symmetry position (using new permutation state)
   *
   *   and then compare dihedral values between the composites. This scheme
   *   has works even if the fused position has changed.
   * - Since newComposite may have reordered the OrientationStates, we have
   *   to be careful which part of the DihedralTuple's we extract symmetry
   *   positions from.
   * - Inversions of the dihedral defining sequence do not invert the sign of
   *   the dihedral
   */
  bool modifiedOrientationIsFirstInNewComposite = (
    newComposite.orientations().first.identifier
    == possiblyModifiedOrientation.identifier
  );

  using DihedralTuple = stereopermutation::Composite::DihedralTuple;
  // This permutator is assigned since that is ensured a few lines earlier
  const std::vector<DihedralTuple>& oldDihedralList = _composite.dihedrals(
    _assignment.value()
  );

  auto getNewSymmetryPosition = [&](unsigned oldSymmetryPosition) -> unsigned {
    const unsigned oldLigandIndex = SymmetryMapHelper::getLigandIndexAt(
      oldSymmetryPosition,
      oldPermutationState.symmetryPositionMap
    );

    const std::vector<AtomIndex>& oldLigand = oldRanking.ligands.at(oldLigandIndex);

    // We assume here that atom indices making up ligands are sorted
    assert(std::is_sorted(std::begin(oldLigand), std::end(oldLigand)));

    /* Find this ligand in the new ranking
     * - In case there is truly only a ranking change (not a rearrangement in
     *   terms of haptic ligands or so), then there should be a vector with
     *   exactly the same elements.
     */
    const auto& newRankingLigands = newPermutator.getRanking().ligands;
    const auto findLigandIter = std::find(
      std::begin(newRankingLigands),
      std::end(newRankingLigands),
      oldLigand
    );

    // TODO this might throw in some cases
    assert(findLigandIter != std::end(newRankingLigands));

    const unsigned newLigandIndex = findLigandIter - std::begin(newRankingLigands);

    return SymmetryMapHelper::getSymmetryPositionOf(
      newLigandIndex,
      newPermutator.getSymmetryPositionMap()
    );
  };

  // Map the set of old dihedral tuples into the new space
  std::vector<DihedralTuple> newCompositeDihedrals;
  newCompositeDihedrals.reserve(oldDihedralList.size());
  for(const DihedralTuple& oldDihedral : oldDihedralList) {
    const unsigned changedSymmetryPosition = detail::select(
      oldDihedral,
      changedIsFirstInOldOrientations
    );

    const unsigned unchangedSymmetryPosition = detail::select(
      oldDihedral,
      !changedIsFirstInOldOrientations
    );

    const unsigned newSymmetryPosition = getNewSymmetryPosition(changedSymmetryPosition);

    if(modifiedOrientationIsFirstInNewComposite) {
      newCompositeDihedrals.emplace_back(
        newSymmetryPosition,
        unchangedSymmetryPosition,
        std::get<2>(oldDihedral)
      );
    } else {
      newCompositeDihedrals.emplace_back(
        unchangedSymmetryPosition,
        newSymmetryPosition,
        std::get<2>(oldDihedral)
      );
    }
  }

  /* Sort the dihedrals for easy comparison. Composite already generates its
   * dihedrals in a sorted fashion, so newComposite's individual permutations'
   * dihedrals need not be sorted.
   */
  std::sort(
    std::begin(newCompositeDihedrals),
    std::end(newCompositeDihedrals)
  );

  /* Find the matching stereopermutation by finding a bijective mapping from
   * all old dihedral tuples of the old stereopermutation to all dihedral
   * tuples of a stereopermutation within the new Composite. Since all lists of
   * DihedralTuples are sorted, we can use the lexicographical vector
   * comparison operator that in turn calls the lexicographical tuple
   * comparison operator:
   *
   * However, the match may not be floating-point exact, and can have pi
   * periodicities!
   */
  auto matchIter = std::find_if(
    std::begin(newComposite),
    std::end(newComposite),
    [&](const std::vector<DihedralTuple>& dihedrals) -> bool {
      return std::lexicographical_compare(
        std::begin(newCompositeDihedrals),
        std::end(newCompositeDihedrals),
        std::begin(dihedrals),
        std::end(dihedrals),
        [](const DihedralTuple& a, const DihedralTuple& b) -> bool {
          return (
            std::get<0>(a) == std::get<0>(b)
            && std::get<1>(a) == std::get<1>(b)
            && detail::piPeriodicFPCompare(std::get<2>(a), std::get<2>(b))
          );
        }
      );
    }
  );

  // There should always be a match
  if(matchIter == std::end(newComposite)) {
    throw std::logic_error("Bug: no match found in new composite.");
  }

  // Overwrite class state
  _assignment = matchIter - std::begin(newComposite);
  _composite = std::move(newComposite);
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
bool BondStereopermutator::Impl::operator < (const Impl& other) const {
  if(_composite < other._composite) {
    return true;
  }

  return assigned() < other.assigned();
}

bool BondStereopermutator::Impl::operator == (const Impl& other) const {
  return (
    _composite == other._composite
    && assigned() == other.assigned()
  );
}

bool BondStereopermutator::Impl::operator != (const Impl& other) const {
  return !(*this == other);
}


} // namespace molassembler

} // namespace Scine
