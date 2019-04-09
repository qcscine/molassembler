/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/AtomStereopermutatorImpl.h"

namespace Scine {

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

Symmetry::Name AtomStereopermutator::up(const Symmetry::Name symmetryName) {
  return Impl::up(symmetryName);
}

Symmetry::Name AtomStereopermutator::down(const Symmetry::Name symmetryName, const unsigned removedSymmetryPosition) {
  return Impl::down(symmetryName, removedSymmetryPosition);
}

void AtomStereopermutator::assign(boost::optional<unsigned> assignment) {
  _pImpl->assign(std::move(assignment));
}

void AtomStereopermutator::assignRandom() {
  _pImpl->assignRandom();
}

void AtomStereopermutator::applyPermutation(const std::vector<AtomIndex>& permutation) {
  _pImpl->applyPermutation(permutation);
}

void AtomStereopermutator::fit(
  const OuterGraph& graph,
  const AngstromWrapper& angstromWrapper
) {
  _pImpl->fit(graph, angstromWrapper);
}

boost::optional<AtomStereopermutator::PropagatedState> AtomStereopermutator::propagate(
  const OuterGraph& graph,
  RankingInformation newRanking,
  boost::optional<Symmetry::Name> symmetryOption
) {
  return _pImpl->propagate(
    graph,
    std::move(newRanking),
    std::move(symmetryOption)
  );
}

void AtomStereopermutator::propagateVertexRemoval(const AtomIndex removedIndex) {
  _pImpl->propagateVertexRemoval(removedIndex);
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
> AtomStereopermutator::minimalChiralConstraints(bool enforce) const {
  return _pImpl->minimalChiralConstraints(enforce);
}

std::string AtomStereopermutator::info() const {
  return _pImpl->info();
}

std::string AtomStereopermutator::rankInfo() const {
  return _pImpl->rankInfo();
}

const PermutationState& AtomStereopermutator::getPermutationState() const {
  return _pImpl->getPermutationState();
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

} // namespace Scine
