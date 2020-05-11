/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutators/BondStereopermutatorImpl.h"

namespace Scine {
namespace Molassembler {

constexpr double BondStereopermutator::assignmentAcceptanceParameter;

/* Impl */

/* Static Impl member functions */

/* BondStereopermutator implementation */
BondStereopermutator::BondStereopermutator(BondStereopermutator&& other) noexcept = default;
BondStereopermutator& BondStereopermutator::operator = (BondStereopermutator&& other) noexcept = default;

BondStereopermutator::BondStereopermutator(const BondStereopermutator& other) : pImpl_(
  std::make_unique<Impl>(*other.pImpl_)
) {}
BondStereopermutator& BondStereopermutator::operator = (const BondStereopermutator& other) {
  *pImpl_ = *other.pImpl_;
  return *this;
}
BondStereopermutator::~BondStereopermutator() = default;

BondStereopermutator::BondStereopermutator(
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const BondIndex& edge,
  const Alignment alignment
) {
  pImpl_ = std::make_unique<Impl>(
    stereopermutatorA,
    stereopermutatorB,
    edge,
    alignment
  );
}

BondStereopermutator::BondStereopermutator(
  const PrivateGraph& graph,
  const StereopermutatorList& stereopermutators,
  const BondIndex& edge,
  const Alignment alignment
) {
  pImpl_ = std::make_unique<Impl>(
    graph,
    stereopermutators,
    edge,
    alignment
  );
}

void BondStereopermutator::assign(boost::optional<unsigned> assignment) {
  pImpl_->assign(std::move(assignment));
}

void BondStereopermutator::assignRandom(Random::Engine& engine) {
  pImpl_->assignRandom(engine);
}

void BondStereopermutator::applyPermutation(const std::vector<AtomIndex>& permutation) {
  pImpl_->applyPermutation(permutation);
}

void BondStereopermutator::fit(
  const AngstromPositions& angstromWrapper,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const FittingMode mode
) {
  pImpl_->fit(
    angstromWrapper,
    stereopermutatorA,
    stereopermutatorB,
    mode
  );
}

void BondStereopermutator::propagateGraphChange(
  const AtomStereopermutatorPropagatedState& oldPermutator,
  const AtomStereopermutator& newPermutator,
  const PrivateGraph& inner,
  const StereopermutatorList& permutators
) {
  pImpl_->propagateGraphChange(
    oldPermutator,
    newPermutator,
    inner,
    permutators
  );
}

BondStereopermutator::Alignment BondStereopermutator::alignment() const {
  return pImpl_->alignment();
}

boost::optional<unsigned> BondStereopermutator::assigned() const {
  return pImpl_->assigned();
}

const Stereopermutations::Composite& BondStereopermutator::composite() const {
  return pImpl_->composite();
}

double BondStereopermutator::dihedral(
  const AtomStereopermutator& stereopermutatorA,
  const SiteIndex siteIndexA,
  const AtomStereopermutator& stereopermutatorB,
  const SiteIndex siteIndexB
) const {
  return pImpl_->dihedral(
    stereopermutatorA,
    siteIndexA,
    stereopermutatorB,
    siteIndexB
  );
}

bool BondStereopermutator::hasSameCompositeOrientation(const BondStereopermutator& other) const {
  return pImpl_->hasSameCompositeOrientation(*other.pImpl_);
}

boost::optional<unsigned> BondStereopermutator::indexOfPermutation() const {
  return pImpl_->indexOfPermutation();
}

unsigned BondStereopermutator::numAssignments() const {
  return pImpl_->numAssignments();
}

unsigned BondStereopermutator::numStereopermutations() const {
  return pImpl_->numStereopermutations();
}

std::string BondStereopermutator::info() const {
  return pImpl_->info();
}

std::string BondStereopermutator::rankInfo() const {
  return pImpl_->rankInfo();
}

BondIndex BondStereopermutator::placement() const {
  return pImpl_->placement();
}

bool BondStereopermutator::operator < (const BondStereopermutator& other) const {
  return *pImpl_ < *other.pImpl_;
}

bool BondStereopermutator::operator == (const BondStereopermutator& other) const {
  return *pImpl_ == *other.pImpl_;
}

bool BondStereopermutator::operator != (const BondStereopermutator& other) const {
  return !(
    *pImpl_ == *other.pImpl_
  );
}

} // namespace Molassembler
} // namespace Scine
