/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "DirectedConformerGenerator.h"
#include "DistanceGeometry/DirectedConformerGeneratorImpl.h"

#include "molassembler/Molecule.h"
#include "molassembler/BondStereopermutator.h"

#include "boost/variant.hpp"

namespace Scine {
namespace molassembler {

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
) : _pImpl(std::make_unique<Impl>(std::move(molecule), bondsToConsider)) {}

DirectedConformerGenerator::DirectedConformerGenerator(DirectedConformerGenerator&& other) noexcept = default;
DirectedConformerGenerator& DirectedConformerGenerator::operator = (DirectedConformerGenerator&& other) noexcept = default;

DirectedConformerGenerator::~DirectedConformerGenerator() = default;

DirectedConformerGenerator::DecisionList DirectedConformerGenerator::generateNewDecisionList() {
  return _pImpl->generateNewDecisionList();
}

bool DirectedConformerGenerator::insert(const DecisionList& decisionList) {
  return _pImpl->insert(decisionList);
}

bool DirectedConformerGenerator::contains(const DecisionList& decisionList) {
  return _pImpl->contains(decisionList);
}

const DirectedConformerGenerator::BondList& DirectedConformerGenerator::bondList() const {
  return _pImpl->bondList();
}

unsigned DirectedConformerGenerator::decisionListSetSize() const {
  return _pImpl->decisionListSetSize();
}

unsigned DirectedConformerGenerator::idealEnsembleSize() const {
  return _pImpl->idealEnsembleSize();
}

outcome::result<Utils::PositionCollection>
DirectedConformerGenerator::generateConformation(
  const DecisionList& decisionList,
  const DistanceGeometry::Configuration& configuration
) {
  return _pImpl->generateConformation(decisionList, configuration);
}

const Molecule& DirectedConformerGenerator::conformationMolecule(const DecisionList& decisionList) {
  return _pImpl->conformationMolecule(decisionList);
}

DirectedConformerGenerator::DecisionList DirectedConformerGenerator::getDecisionList(Utils::PositionCollection positions) const {
  return _pImpl->getDecisionList(std::move(positions));
}

} // namespace molassembler
} // namespace Scine