// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Stereocenters/AtomStereocenterImpl.h"

namespace molassembler {

/* AtomStereocenter implementations */
AtomStereocenter::AtomStereocenter(
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

AtomStereocenter::AtomStereocenter(AtomStereocenter&& other) noexcept = default;
AtomStereocenter& AtomStereocenter::operator = (AtomStereocenter&& other) noexcept = default;

AtomStereocenter::AtomStereocenter(const AtomStereocenter& other) : _pImpl(
  std::make_unique<Impl>(*other._pImpl)
) {}
AtomStereocenter& AtomStereocenter::operator = (const AtomStereocenter& other) {
  *_pImpl = *other._pImpl;
  return *this;
}

AtomStereocenter::~AtomStereocenter() = default;

void AtomStereocenter::addSubstituent(
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

void AtomStereocenter::assign(boost::optional<unsigned> assignment) {
  _pImpl->assign(std::move(assignment));
}

void AtomStereocenter::assignRandom() {
  _pImpl->assignRandom();
}

void AtomStereocenter::fit(
  const OuterGraph& graph,
  const AngstromWrapper& angstromWrapper,
  const std::vector<Symmetry::Name>& excludeSymmetries
) {
  _pImpl->fit(graph, angstromWrapper, excludeSymmetries);
}

void AtomStereocenter::propagateGraphChange(
  const OuterGraph& graph,
  RankingInformation newRanking
) {
  _pImpl->propagateGraphChange(
    graph,
    std::move(newRanking)
  );
}

void AtomStereocenter::propagateVertexRemoval(const AtomIndex removedIndex) {
  _pImpl->propagateVertexRemoval(removedIndex);
}

void AtomStereocenter::removeSubstituent(
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

void AtomStereocenter::setSymmetry(
  const Symmetry::Name symmetryName,
  const OuterGraph& graph
) {
  _pImpl->setSymmetry(symmetryName, graph);
}

/* Information */
double AtomStereocenter::angle(
  const unsigned i,
  const unsigned j
) const {
  return _pImpl->angle(i, j);
}

boost::optional<unsigned> AtomStereocenter::assigned() const {
  return _pImpl->assigned();
}

AtomIndex AtomStereocenter::centralIndex() const {
  return _pImpl->centralIndex();
}

boost::optional<unsigned> AtomStereocenter::indexOfPermutation() const {
  return _pImpl->indexOfPermutation();
}

std::vector<
  std::array<boost::optional<unsigned>, 4>
> AtomStereocenter::minimalChiralityConstraints() const {
  return _pImpl->minimalChiralityConstraints();
}

std::vector<DistanceGeometry::ChiralityConstraint> AtomStereocenter::chiralityConstraints(
  double looseningMultiplier
) const {
  return _pImpl->chiralityConstraints(looseningMultiplier);
}

std::string AtomStereocenter::info() const {
  return _pImpl->info();
}

std::string AtomStereocenter::rankInfo() const {
  return _pImpl->rankInfo();
}

const RankingInformation& AtomStereocenter::getRanking() const {
  return _pImpl->getRanking();
}

Symmetry::Name AtomStereocenter::getSymmetry() const {
  return _pImpl->getSymmetry();
}

std::vector<unsigned> AtomStereocenter::getSymmetryPositionMap() const {
  return _pImpl->getSymmetryPositionMap();
}

unsigned AtomStereocenter::numAssignments() const {
  return _pImpl->numAssignments();
}

unsigned AtomStereocenter::numStereopermutations() const {
  return _pImpl->numStereopermutations();
}

void AtomStereocenter::setModelInformation(
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
bool AtomStereocenter::operator == (const AtomStereocenter& other) const {
  return *_pImpl == *other._pImpl;
}

bool AtomStereocenter::operator != (const AtomStereocenter& other) const {
  return !(*_pImpl == *other._pImpl);
}
bool AtomStereocenter::operator < (const AtomStereocenter& other) const {
  return *_pImpl < *other._pImpl;
}

} // namespace molassembler
