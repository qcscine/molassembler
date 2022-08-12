/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "DirectedConformerGenerator.h"
#include "DistanceGeometry/DirectedConformerGeneratorImpl.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/Stereopermutation/Composites.h"
#include "Molassembler/Detail/Cartesian.h"
#include "Molassembler/Temple/Adaptors/CyclicFrame.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functor.h"

#include "boost/variant.hpp"

namespace Scine {
namespace Molassembler {

constexpr std::uint8_t DirectedConformerGenerator::unknownDecision;

/* Static functions */
unsigned DirectedConformerGenerator::distance(
  const DecisionList& a,
  const DecisionList& b,
  const DecisionList& bounds
) {
  return Impl::distance(a, b, bounds);
}

boost::variant<DirectedConformerGenerator::IgnoreReason, BondStereopermutator>
DirectedConformerGenerator::considerBond(
  const BondIndex& bondIndex,
  const Molecule& molecule,
  const BondStereopermutator::Alignment alignment
) {
  return Impl::considerBond(bondIndex, molecule, alignment);
}

DirectedConformerGenerator::DirectedConformerGenerator(
  Molecule molecule,
  const BondStereopermutator::Alignment alignment,
  const BondList& bondsToConsider
) : pImpl_(std::make_unique<Impl>(std::move(molecule), alignment, bondsToConsider)) {}

DirectedConformerGenerator::DirectedConformerGenerator(DirectedConformerGenerator&& other) noexcept = default;
DirectedConformerGenerator& DirectedConformerGenerator::operator = (DirectedConformerGenerator&& other) noexcept = default;

DirectedConformerGenerator::~DirectedConformerGenerator() = default;

DirectedConformerGenerator::DecisionList DirectedConformerGenerator::generateNewDecisionList(Random::Engine& engine) {
  return pImpl_->generateNewDecisionList(engine);
}

bool DirectedConformerGenerator::insert(const DecisionList& decisionList) {
  return pImpl_->insert(decisionList);
}

bool DirectedConformerGenerator::contains(const DecisionList& decisionList) const {
  return pImpl_->contains(decisionList);
}

BondStereopermutator::Alignment DirectedConformerGenerator::alignment() const {
  return pImpl_->alignment();
}

const DirectedConformerGenerator::BondList& DirectedConformerGenerator::bondList() const {
  return pImpl_->bondList();
}

unsigned DirectedConformerGenerator::decisionListSetSize() const {
  return pImpl_->decisionListSetSize();
}

unsigned DirectedConformerGenerator::idealEnsembleSize() const {
  return pImpl_->idealEnsembleSize();
}

outcome::result<Utils::PositionCollection>
DirectedConformerGenerator::generateRandomConformation(
  const DecisionList& decisionList,
  const DistanceGeometry::Configuration& configuration,
  const BondStereopermutator::FittingMode fitting
) const {
  return pImpl_->generateRandomConformation(decisionList, configuration, fitting);
}

outcome::result<Utils::PositionCollection>
DirectedConformerGenerator::generateConformation(
  const DecisionList& decisionList,
  const unsigned seed,
  const DistanceGeometry::Configuration& configuration,
  const BondStereopermutator::FittingMode fitting
) const {
  return pImpl_->generateConformation(decisionList, seed, configuration, fitting);
}

Molecule DirectedConformerGenerator::conformationMolecule(const DecisionList& decisionList) const {
  return pImpl_->conformationMolecule(decisionList);
}

DirectedConformerGenerator::DecisionList DirectedConformerGenerator::getDecisionList(
  const Utils::AtomCollection& atomCollection,
  const BondStereopermutator::FittingMode mode
) const {
  return pImpl_->getDecisionList(atomCollection, mode);
}

DirectedConformerGenerator::DecisionList DirectedConformerGenerator::getDecisionList(
  const Utils::PositionCollection& positions,
  const BondStereopermutator::FittingMode mode
) const {
  return pImpl_->getDecisionList(positions, mode);
}

void DirectedConformerGenerator::enumerate(
  std::function<void(const DecisionList&, Utils::PositionCollection)> callback,
  unsigned seed,
  const EnumerationSettings& settings
) {
  return pImpl_->enumerate(std::move(callback), seed, settings);
}

void DirectedConformerGenerator::enumerateRandom(
  std::function<void(const DecisionList&, Utils::PositionCollection)> callback,
  const EnumerationSettings& settings
) {
  Random::Engine& engine = randomnessEngine();
  return pImpl_->enumerate(std::move(callback), engine(), settings);
}

DirectedConformerGenerator::Relabeler DirectedConformerGenerator::relabeler() const {
  return pImpl_->relabeler();
}

std::vector<int>
DirectedConformerGenerator::binMidpointIntegers(const DecisionList& decisions) const {
  return pImpl_->binMidpointIntegers(decisions);
}

std::vector<std::pair<int, int>>
DirectedConformerGenerator::binBounds(const DecisionList& decisions) const {
  return pImpl_->binBounds(decisions);
}

std::pair<double, double> DirectedConformerGenerator::Relabeler::makeBounds(
  double phi,
  double tolerance
) {
  return Temple::map(
    std::make_pair(phi - tolerance, phi + tolerance),
    &Cartesian::signedDihedralAngle
  );
}

std::pair<int, int> DirectedConformerGenerator::Relabeler::integerBounds(
  const std::pair<double, double>& bounds
) {
  return std::make_pair(
    static_cast<int>(std::floor(Temple::Math::toDegrees(bounds.first))),
    static_cast<int>(std::ceil(Temple::Math::toDegrees(bounds.second)))
  );
}

DirectedConformerGenerator::Relabeler::Relabeler(
  const DirectedConformerGenerator::BondList& bonds,
  const Molecule& mol
) {
  // Determine the dominant dihedral sequences at each bond
  for(const BondIndex& bond : bonds) {
    const auto& stereopermutator = mol.stereopermutators().at(bond);
    const auto& composite = stereopermutator.composite();

    const AtomIndex leftPlacement = composite.orientations().first.identifier;
    const AtomIndex rightPlacement = composite.orientations().second.identifier;

    const auto& left = mol.stereopermutators().at(leftPlacement);
    const auto& right = mol.stereopermutators().at(rightPlacement);

    mol.rankPriority(leftPlacement);
    mol.rankPriority(rightPlacement);


    const auto& dihedrals = composite.allPermutations().at(0).dihedrals;
    SiteIndex dominantLeft;
    SiteIndex dominantRight;
    unsigned maxLeftRanking = std::numeric_limits<unsigned>::max();
    unsigned maxRightRanking = std::numeric_limits<unsigned>::max();
    for(const auto& dihedral : dihedrals) {
      const SiteIndex leftSite = left.getShapePositionMap().indexOf(std::get<0>(dihedral));
      const SiteIndex rightSite = right.getShapePositionMap().indexOf(std::get<1>(dihedral));
      const auto leftRanking = left.getRanking().getRankedIndexOfSite(leftSite);
      const auto rightRanking = right.getRanking().getRankedIndexOfSite(rightSite);

      const bool betterLeft = (leftRanking < maxLeftRanking && rightRanking <= maxRightRanking);
      const bool betterRight = (leftRanking <= maxLeftRanking && rightRanking < maxRightRanking);
      if(betterLeft || betterRight) {
        dominantLeft = leftSite;
        dominantRight = rightSite;
        maxLeftRanking = leftRanking;
        maxRightRanking = rightRanking;
      }
    }

    sequences.push_back(
      DihedralInfo {
        left.getRanking().sites.at(dominantLeft),
        leftPlacement,
        rightPlacement,
        right.getRanking().sites.at(dominantRight),
        composite.rotationalAxisSymmetryOrder()
      }
    );
  }

  observedDihedrals.resize(sequences.size());
}

std::vector<double> DirectedConformerGenerator::Relabeler::add(
  const Utils::PositionCollection& positions
) {
  auto structureDihedrals = Temple::map(
    sequences,
    [&](const auto& sequence) -> double {
      double dihedral = Cartesian::dihedral(
        Cartesian::averagePosition(positions, sequence.is),
        positions.row(sequence.j),
        positions.row(sequence.k),
        Cartesian::averagePosition(positions, sequence.ls)
      );

      // Reduce by rotational symmetry if present
      if(sequence.symmetryOrder > 1) {
        const double tmp = 2 * M_PI / sequence.symmetryOrder;
        const double radians = std::fmod(
          Cartesian::positiveDihedralAngle(dihedral),
          tmp
        );
        dihedral = radians - std::floor((radians + tmp / 2.0) / tmp) * tmp;
      }

      return dihedral;
    }
  );

  // Distribute the dihedrals into each sequence's list of dihedrals
  auto dihedralIter = std::begin(observedDihedrals);
  assert(observedDihedrals.size() == structureDihedrals.size());
  for(double dihedral : structureDihedrals) {
    dihedralIter->push_back(dihedral);
    ++dihedralIter;
  }

  return structureDihedrals;
}

DirectedConformerGenerator::Relabeler::Intervals
DirectedConformerGenerator::Relabeler::densityBins(
  const std::vector<double>& dihedrals,
  const double delta,
  unsigned symmetryOrder
) {
  if(dihedrals.empty()) {
    throw std::logic_error("Cannot make bins for empty list of dihedrals");
  }

  if(dihedrals.size() == 1) {
    return Intervals(1, Interval(dihedrals.front(), dihedrals.front()));
  }

  const auto sortedDihedrals = Temple::sorted(dihedrals);
  const double boundary = 2 * M_PI / symmetryOrder;

  std::pair<double, double> interval {
    sortedDihedrals.front(),
    sortedDihedrals.front()
  };
  bool closeInterval = false;

  Intervals intervals;

  Temple::forEach(
    Temple::Adaptors::cyclicFrame<2>(sortedDihedrals),
    [&](const double prev, const double next) {
      if(prev <= next) {
        closeInterval = (next - prev) > delta;
      } else {
        closeInterval = (next + boundary - prev) > delta;
      }

      if(closeInterval) {
        interval.second = prev;
        intervals.push_back(interval);
        interval.first = next;
      }
    }
  );

  if(!closeInterval) {
    if(intervals.empty()) {
      interval.second = sortedDihedrals.back();
      intervals.push_back(interval);
    } else {
      /* Merge the first interval with the remaining open one whose start is
       * stored in interval
       */
      intervals.front().first = interval.first;
    }
  }

  return intervals;
}

std::vector<DirectedConformerGenerator::Relabeler::Intervals>
DirectedConformerGenerator::Relabeler::bins(const double delta) const {
  return Temple::map(
    Temple::Adaptors::zip(observedDihedrals, sequences),
    [&](const auto& dihedrals, const DihedralInfo& sequence) -> Intervals {
      return densityBins(dihedrals, delta, sequence.symmetryOrder);
    }
  );
}

std::vector<std::vector<unsigned>>
DirectedConformerGenerator::Relabeler::binIndices(
  const std::vector<Intervals>& allBins
) const {
  const unsigned structureCount = observedDihedrals.front().size();
  const unsigned bondCount = observedDihedrals.size();

  std::vector<std::vector<unsigned>> relabeling(structureCount, std::vector<unsigned>(bondCount));

#pragma omp parallel for collapse(2)
  for(unsigned structure = 0; structure < structureCount; ++structure) {
    for(unsigned bond = 0; bond < bondCount; ++bond) {
      const auto& bins = allBins.at(bond);
      double observedDihedral = observedDihedrals.at(bond).at(structure);

      const auto findIter = Temple::find_if(
        bins,
        [&](const Interval& interval) -> bool {
          if(interval.first <= interval.second) {
            return interval.first <= observedDihedral && observedDihedral <= interval.second;
          }

          return interval.first <= observedDihedral || observedDihedral <= interval.second;
        }
      );

      assert(findIter != std::end(bins));

      relabeling.at(structure).at(bond) = findIter - std::begin(bins);
    }
  }


  return relabeling;
}

std::vector<std::vector<int>>
DirectedConformerGenerator::Relabeler::binMidpointIntegers(
  const std::vector<std::vector<unsigned>>& binIndices,
  const std::vector<Intervals>& allBins
) const {
  const auto intervalMidpoint = [](
    const Interval& interval,
    const unsigned symmetryOrder
  ) -> int {
    if(interval.first <= interval.second) {
      return static_cast<int>(std::round(180 * (interval.first + interval.second) / (2 * M_PI)));
    }

    const double boundary = 2 * M_PI / symmetryOrder;
    const double average = Cartesian::signedDihedralAngle(
      (interval.first + interval.second + boundary) / 2
    );
    return static_cast<int>(std::round(180 * average / M_PI));
  };

  const auto binMidpointIntegers = Temple::map(
    Temple::Adaptors::zip(allBins, sequences),
    [&](const auto& intervals, const DihedralInfo& sequence) {
      return Temple::map(
        intervals,
        [&](const Interval& v) {
          return intervalMidpoint(v, sequence.symmetryOrder);
        }
      );
    }
  );

  const unsigned structureCount = binIndices.size();
  const unsigned bondCount = binIndices.front().size();

  std::vector<std::vector<int>> relabeling(structureCount, std::vector<int>(bondCount));

#pragma omp parallel for collapse(2)
  for(unsigned structure = 0; structure < structureCount; ++structure) {
    for(unsigned bond = 0; bond < bondCount; ++bond) {
      relabeling.at(structure).at(bond) = binMidpointIntegers.at(bond).at(
        binIndices.at(structure).at(bond)
      );
    }
  }

  return relabeling;
}

std::vector<std::vector<std::pair<int, int>>>
DirectedConformerGenerator::Relabeler::binBounds(
  const std::vector<std::vector<unsigned>>& binIndices,
  const std::vector<Intervals>& allBins
) const {
  const auto binBoundIntegers = Temple::map(
    Temple::Adaptors::zip(allBins, sequences),
    [&](const auto& intervals, const DihedralInfo& sequence) {
      return Temple::map(
        intervals,
        [&](const Interval& v) -> std::pair<int, int> {
          const auto reducedBounds = Temple::map(
            v,
            [&](const double phi) -> double {
              // Reduce the dihedral angle by the symmetry order
              return Cartesian::signedDihedralAngle(
                std::fmod(
                  Cartesian::positiveDihedralAngle(phi),
                  2 * M_PI / sequence.symmetryOrder
                )
              );
            }
          );

          return std::make_pair(
            static_cast<int>(std::floor(Temple::Math::toDegrees(reducedBounds.first))),
            static_cast<int>(std::ceil(Temple::Math::toDegrees(reducedBounds.second)))
          );
        }
      );
    }
  );

  const unsigned structureCount = binIndices.size();
  const unsigned bondCount = binIndices.front().size();

  std::vector<std::vector<std::pair<int, int>>> relabeling(structureCount, std::vector<std::pair<int, int>>(bondCount));

#pragma omp parallel for collapse(2)
  for(unsigned structure = 0; structure < structureCount; ++structure) {
    for(unsigned bond = 0; bond < bondCount; ++bond) {
      relabeling.at(structure).at(bond) = binBoundIntegers.at(bond).at(
        binIndices.at(structure).at(bond)
      );
    }
  }

  return relabeling;
}

} // namespace Molassembler
} // namespace Scine
