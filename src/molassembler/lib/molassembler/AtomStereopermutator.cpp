/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/AtomStereopermutatorImpl.h"

namespace Scine {

namespace molassembler {

/* AtomStereopermutator implementations */
AtomStereopermutator::AtomStereopermutator(
  const Graph& graph,
  const shapes::Shape shape,
  const AtomIndex centerAtom,
  RankingInformation ranking
) : _pImpl(
  std::make_unique<Impl>(
    graph,
    shape,
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

shapes::Shape AtomStereopermutator::up(const shapes::Shape shape) {
  return Impl::up(shape);
}

shapes::Shape AtomStereopermutator::down(const shapes::Shape shape, const unsigned removedShapePosition) {
  return Impl::down(shape, removedShapePosition);
}

void AtomStereopermutator::assign(boost::optional<unsigned> assignment) {
  _pImpl->assign(std::move(assignment));
}

void AtomStereopermutator::assignRandom(random::Engine& engine) {
  _pImpl->assignRandom(engine);
}

void AtomStereopermutator::applyPermutation(const std::vector<AtomIndex>& permutation) {
  _pImpl->applyPermutation(permutation);
}

void AtomStereopermutator::fit(
  const Graph& graph,
  const AngstromWrapper& angstromWrapper
) {
  _pImpl->fit(graph, angstromWrapper);
}

boost::optional<AtomStereopermutator::PropagatedState> AtomStereopermutator::propagate(
  const Graph& graph,
  RankingInformation newRanking,
  boost::optional<shapes::Shape> shapeOption
) {
  return _pImpl->propagate(
    graph,
    std::move(newRanking),
    std::move(shapeOption)
  );
}

void AtomStereopermutator::propagateVertexRemoval(const AtomIndex removedIndex) {
  _pImpl->propagateVertexRemoval(removedIndex);
}

void AtomStereopermutator::setShape(
  const shapes::Shape shape,
  const Graph& graph
) {
  _pImpl->setShape(shape, graph);
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

std::vector<AtomStereopermutator::MinimalChiralConstraint>
AtomStereopermutator::minimalChiralConstraints(bool enforce) const {
  return _pImpl->minimalChiralConstraints(enforce);
}

std::string AtomStereopermutator::info() const {
  return _pImpl->info();
}

std::string AtomStereopermutator::rankInfo() const {
  return _pImpl->rankInfo();
}

const stereopermutators::Abstract& AtomStereopermutator::getAbstract() const {
  return _pImpl->getAbstract();
}

const stereopermutators::Feasible& AtomStereopermutator::getFeasible() const {
  return _pImpl->getFeasible();
}

const RankingInformation& AtomStereopermutator::getRanking() const {
  return _pImpl->getRanking();
}

shapes::Shape AtomStereopermutator::getShape() const {
  return _pImpl->getShape();
}

const std::vector<unsigned>& AtomStereopermutator::getShapePositionMap() const {
  return _pImpl->getShapePositionMap();
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
