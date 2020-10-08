/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
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
  return pImpl_->enumerate(callback, seed, settings);
}

void DirectedConformerGenerator::enumerateRandom(
  std::function<void(const DecisionList&, Utils::PositionCollection)> callback,
  const EnumerationSettings& settings
) {
  Random::Engine& engine = randomnessEngine();
  return pImpl_->enumerate(callback, engine(), settings);
}

DirectedConformerGenerator::Relabeler DirectedConformerGenerator::relabeler() const {
  return pImpl_->relabeler();
}

std::vector<int>
DirectedConformerGenerator::binMidpointIntegers(
  const DecisionList& decisions
) const {
  return pImpl_->binMidpointIntegers(decisions);
}

DirectedConformerGenerator::Relabeler::Relabeler(
  DirectedConformerGenerator::BondList bonds,
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

    const auto& dominantDihedralTuple = composite.allPermutations().at(0).dihedrals.front();
    const SiteIndex leftSite = left.getShapePositionMap().indexOf(std::get<0>(dominantDihedralTuple));
    const SiteIndex rightSite = right.getShapePositionMap().indexOf(std::get<1>(dominantDihedralTuple));

    sequences.push_back(
      DihedralInfo {
        left.getRanking().sites.at(leftSite),
        leftPlacement,
        rightPlacement,
        right.getRanking().sites.at(rightSite),
        composite.rotationalAxisSymmetryOrder()
      }
    );
  }

  observedDihedrals.resize(sequences.size());
}

void DirectedConformerGenerator::Relabeler::add(
  const Utils::PositionCollection& positions
) {
  auto sequenceIter = std::cbegin(sequences);
  const auto bondEnd = std::cend(sequences);
  auto dihedralIter = std::begin(observedDihedrals);
  const auto dihedralEnd = std::end(observedDihedrals);

  while(sequenceIter != bondEnd && dihedralIter != dihedralEnd) {
    double dihedral = Cartesian::dihedral(
      Cartesian::averagePosition(positions, sequenceIter->is),
      positions.row(sequenceIter->j),
      positions.row(sequenceIter->k),
      Cartesian::averagePosition(positions, sequenceIter->ls)
    );

    // Reduce by rotational symmetry if present
    if(sequenceIter->symmetryOrder > 1) {
      dihedral = Cartesian::signedDihedralAngle(
        std::fmod(
          Cartesian::positiveDihedralAngle(dihedral),
          2 * M_PI / sequenceIter->symmetryOrder
        )
      );
    }

    dihedralIter->push_back(dihedral);

    ++sequenceIter;
    ++dihedralIter;
  }
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
  bool closeInterval;

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
    unsigned symmetryOrder
  ) -> int {
    if(interval.first <= interval.second) {
      return std::round(180 * (interval.first + interval.second) / (2 * M_PI));
    }

    double boundary = 2 * M_PI / symmetryOrder;
    double average = Cartesian::signedDihedralAngle(
      (interval.first + interval.second + boundary) / 2
    );
    return std::round(180 * average / M_PI);
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

} // namespace Molassembler
} // namespace Scine
