// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Stereopermutators/BondStereopermutatorImpl.h"

namespace molassembler {

constexpr double BondStereopermutator::chiralityConstraintTolerance;
constexpr double BondStereopermutator::assignmentAcceptanceDihedralThreshold;

/* Impl */

/* Static Impl member functions */

/* BondStereopermutator implementation */
BondStereopermutator::BondStereopermutator(BondStereopermutator&& other) noexcept = default;
BondStereopermutator& BondStereopermutator::operator = (BondStereopermutator&& other) noexcept = default;

BondStereopermutator::BondStereopermutator(const BondStereopermutator& other) : _pImpl(
  std::make_unique<Impl>(*other._pImpl)
) {}
BondStereopermutator& BondStereopermutator::operator = (const BondStereopermutator& other) {
  *_pImpl = *other._pImpl;
  return *this;
}
BondStereopermutator::~BondStereopermutator() = default;

BondStereopermutator::BondStereopermutator(
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB,
  const BondIndex& edge
) {
  _pImpl = std::make_unique<Impl>(
    stereopermutatorA,
    stereopermutatorB,
    edge
  );
}

void BondStereopermutator::assign(boost::optional<unsigned> assignment) {
  _pImpl -> assign(std::move(assignment));
}

void BondStereopermutator::assignRandom() {
  _pImpl -> assignRandom();
}

void BondStereopermutator::fit(
  const AngstromWrapper& angstromWrapper,
  const AtomStereopermutator& stereopermutatorA,
  const AtomStereopermutator& stereopermutatorB
) {
  _pImpl -> fit(
    angstromWrapper,
    stereopermutatorA,
    stereopermutatorB
  );
}

void BondStereopermutator::propagateGraphChange(
  const AtomStereopermutator& oldPermutator,
  const AtomStereopermutator& newPermutator
) {
  _pImpl -> propagateGraphChange(oldPermutator, newPermutator);
}

boost::optional<unsigned> BondStereopermutator::assigned() const {
  return _pImpl -> assigned();
}

const stereopermutation::Composite& BondStereopermutator::composite() const {
  return _pImpl -> composite();
}

bool BondStereopermutator::hasSameCompositeOrientation(const BondStereopermutator& other) const {
  return _pImpl -> hasSameCompositeOrientation(*other._pImpl);
}

boost::optional<unsigned> BondStereopermutator::indexOfPermutation() const {
  return _pImpl -> indexOfPermutation();
}

unsigned BondStereopermutator::numAssignments() const {
  return _pImpl -> numAssignments();
}

unsigned BondStereopermutator::numStereopermutations() const {
  return _pImpl -> numStereopermutations();
}

std::string BondStereopermutator::info() const {
  return _pImpl -> info();
}

std::string BondStereopermutator::rankInfo() const {
  return _pImpl -> rankInfo();
}

BondIndex BondStereopermutator::edge() const {
  return _pImpl -> edge();
}

bool BondStereopermutator::operator == (const BondStereopermutator& other) const {
  return *_pImpl == *other._pImpl;
}

bool BondStereopermutator::operator != (const BondStereopermutator& other) const {
  return !(
    *_pImpl == *other._pImpl
  );
}

} // namespace molassembler
