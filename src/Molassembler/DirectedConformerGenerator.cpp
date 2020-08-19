/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "DirectedConformerGenerator.h"
#include "DistanceGeometry/DirectedConformerGeneratorImpl.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/BondStereopermutator.h"

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
  const std::unordered_map<AtomIndex, unsigned>& smallestCycleMap
) {
  return Impl::considerBond(
    bondIndex,
    molecule,
    smallestCycleMap
  );
}

DirectedConformerGenerator::DirectedConformerGenerator(
  Molecule molecule,
  const BondList& bondsToConsider
) : pImpl_(std::make_unique<Impl>(std::move(molecule), bondsToConsider)) {}

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

} // namespace Molassembler
} // namespace Scine
