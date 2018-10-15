// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Stereopermutators/AtomStereopermutatorImpl.h"

namespace molassembler {

/* AtomStereopermutator implementations */
AtomStereopermutator::AtomStereopermutator(
  const OuterGraph& graph,
  const Symmetry::Name symmetry,
  const AtomIndex centerAtom,
  RankingInformation ranking
) : _pImpl(
  std::make_unique<Impl>(
    graph,
    symmetry,
    centerAtom,
    std::move(ranking)
  )
) {}

AtomStereopermutator::AtomStereopermutator(AtomStereopermutator&& other) noexcept = default;
AtomStereopermutator& AtomStereopermutator::operator = (AtomStereopermutator&& other) noexcept = default;

AtomStereopermutator::AtomStereopermutator(const AtomStereopermutator& other) : _pImpl(
  std::make_unique<Impl>(*other._pImpl)
) {}
AtomStereopermutator& AtomStereopermutator::operator = (const AtomStereopermutator& other) {
  *_pImpl = *other._pImpl;
  return *this;
}

AtomStereopermutator::~AtomStereopermutator() = default;

void AtomStereopermutator::addSubstituent(
  const OuterGraph& graph,
  const AtomIndex newSubstituentIndex,
  RankingInformation newRanking,
  const Symmetry::Name newSymmetry,
  const ChiralStatePreservation preservationOption
) {
  _pImpl->addSubstituent(
    graph,
    newSubstituentIndex,
    std::move(newRanking),
    newSymmetry,
    preservationOption
  );
}

void AtomStereopermutator::assign(boost::optional<unsigned> assignment) {
  _pImpl->assign(std::move(assignment));
}

void AtomStereopermutator::assignRandom() {
  _pImpl->assignRandom();
}

void AtomStereopermutator::fit(
  const OuterGraph& graph,
  const AngstromWrapper& angstromWrapper,
  const std::vector<Symmetry::Name>& excludeSymmetries
) {
  _pImpl->fit(graph, angstromWrapper, excludeSymmetries);
}

void AtomStereopermutator::propagateGraphChange(
  const OuterGraph& graph,
  RankingInformation newRanking
) {
  _pImpl->propagateGraphChange(
    graph,
    std::move(newRanking)
  );
}

void AtomStereopermutator::propagateVertexRemoval(const AtomIndex removedIndex) {
  _pImpl->propagateVertexRemoval(removedIndex);
}

void AtomStereopermutator::removeSubstituent(
  const OuterGraph& graph,
  const AtomIndex which,
  RankingInformation newRanking,
  const Symmetry::Name newSymmetry,
  const ChiralStatePreservation preservationOption
) {
  _pImpl->removeSubstituent(
    graph,
    which,
    std::move(newRanking),
    newSymmetry,
    preservationOption
  );
}

void AtomStereopermutator::setSymmetry(
  const Symmetry::Name symmetryName,
  const OuterGraph& graph
) {
  _pImpl->setSymmetry(symmetryName, graph);
}

/* Information */
double AtomStereopermutator::angle(
  const unsigned i,
  const unsigned j
) const {
  return _pImpl->angle(i, j);
}

boost::optional<unsigned> AtomStereopermutator::assigned() const {
  return _pImpl->assigned();
}

AtomIndex AtomStereopermutator::centralIndex() const {
  return _pImpl->centralIndex();
}

boost::optional<unsigned> AtomStereopermutator::indexOfPermutation() const {
  return _pImpl->indexOfPermutation();
}

std::vector<
  std::array<boost::optional<unsigned>, 4>
> AtomStereopermutator::minimalChiralityConstraints() const {
  return _pImpl->minimalChiralityConstraints();
}

std::vector<DistanceGeometry::ChiralityConstraint> AtomStereopermutator::chiralityConstraints(
  double looseningMultiplier
) const {
  return _pImpl->chiralityConstraints(looseningMultiplier);
}

std::string AtomStereopermutator::info() const {
  return _pImpl->info();
}

std::string AtomStereopermutator::rankInfo() const {
  return _pImpl->rankInfo();
}

const RankingInformation& AtomStereopermutator::getRanking() const {
  return _pImpl->getRanking();
}

Symmetry::Name AtomStereopermutator::getSymmetry() const {
  return _pImpl->getSymmetry();
}

std::vector<unsigned> AtomStereopermutator::getSymmetryPositionMap() const {
  return _pImpl->getSymmetryPositionMap();
}

unsigned AtomStereopermutator::numAssignments() const {
  return _pImpl->numAssignments();
}

unsigned AtomStereopermutator::numStereopermutations() const {
  return _pImpl->numStereopermutations();
}

void AtomStereopermutator::setModelInformation(
  DistanceGeometry::SpatialModel& model,
  const std::function<double(const AtomIndex)>& cycleMultiplierForIndex,
  const double looseningMultiplier
) const {
  _pImpl->setModelInformation(
    model,
    cycleMultiplierForIndex,
    looseningMultiplier
  );
}


/* Operators */
bool AtomStereopermutator::operator == (const AtomStereopermutator& other) const {
  return *_pImpl == *other._pImpl;
}

bool AtomStereopermutator::operator != (const AtomStereopermutator& other) const {
  return !(*_pImpl == *other._pImpl);
}
bool AtomStereopermutator::operator < (const AtomStereopermutator& other) const {
  return *_pImpl < *other._pImpl;
}

} // namespace molassembler
