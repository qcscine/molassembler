/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutators/AtomStereopermutatorImpl.h"

#include "Molassembler/Graph.h"

namespace Scine {
namespace Molassembler {

/* AtomStereopermutator implementations */
AtomStereopermutator::AtomStereopermutator(
  const Graph& graph,
  const Shapes::Shape shape,
  const AtomIndex centerAtom,
  RankingInformation ranking
) : pImpl_(
  std::make_unique<Impl>(
    graph,
    shape,
    centerAtom,
    std::move(ranking)
  )
) {}

AtomStereopermutator::AtomStereopermutator(AtomStereopermutator&& other) noexcept = default;
AtomStereopermutator& AtomStereopermutator::operator = (AtomStereopermutator&& other) noexcept = default;

AtomStereopermutator::AtomStereopermutator(const AtomStereopermutator& other) : pImpl_(
  std::make_unique<Impl>(*other.pImpl_)
) {}
AtomStereopermutator& AtomStereopermutator::operator = (const AtomStereopermutator& other) {
  *pImpl_ = *other.pImpl_;
  return *this;
}

AtomStereopermutator::~AtomStereopermutator() = default;

Shapes::Shape AtomStereopermutator::up(const Shapes::Shape shape) {
  return Impl::up(shape);
}

Shapes::Shape AtomStereopermutator::down(const Shapes::Shape shape, const Shapes::Vertex removedVertex) {
  return Impl::down(shape, removedVertex);
}

void AtomStereopermutator::assign(boost::optional<unsigned> assignment) {
  pImpl_->assign(std::move(assignment));
}

void AtomStereopermutator::assignRandom(Random::Engine& engine) {
  pImpl_->assignRandom(engine);
}

void AtomStereopermutator::applyPermutation(const std::vector<AtomIndex>& permutation) {
  pImpl_->applyPermutation(permutation);
}

boost::optional<AtomStereopermutator::ShapeMap> AtomStereopermutator::fit(
  const Graph& graph,
  const AngstromPositions& angstromWrapper
) {
  return pImpl_->fit(graph, angstromWrapper);
}

boost::optional<AtomStereopermutator::PropagatedState> AtomStereopermutator::propagate(
  const Graph& graph,
  RankingInformation newRanking,
  boost::optional<Shapes::Shape> shapeOption
) {
  return pImpl_->propagate(
    graph,
    std::move(newRanking),
    std::move(shapeOption)
  );
}

void AtomStereopermutator::propagateVertexRemoval(const AtomIndex removedIndex) {
  pImpl_->propagateVertexRemoval(removedIndex);
}

void AtomStereopermutator::setShape(
  const Shapes::Shape shape,
  const Graph& graph
) {
  pImpl_->setShape(shape, graph);
}

void AtomStereopermutator::thermalize(const bool thermalization) {
  pImpl_->thermalize(thermalization);
}

/* Information */
double AtomStereopermutator::angle(
  const SiteIndex i,
  const SiteIndex j
) const {
  return pImpl_->angle(i, j);
}

boost::optional<unsigned> AtomStereopermutator::assigned() const {
  return pImpl_->assigned();
}

AtomIndex AtomStereopermutator::placement() const {
  return pImpl_->placement();
}

boost::optional<unsigned> AtomStereopermutator::indexOfPermutation() const {
  return pImpl_->indexOfPermutation();
}

std::vector<AtomStereopermutator::MinimalChiralConstraint>
AtomStereopermutator::minimalChiralConstraints(bool enforce) const {
  return pImpl_->minimalChiralConstraints(enforce);
}

std::string AtomStereopermutator::info() const {
  return pImpl_->info();
}

std::string AtomStereopermutator::rankInfo() const {
  return pImpl_->rankInfo();
}

std::vector<std::vector<SiteIndex>> AtomStereopermutator::siteGroups() const {
  return pImpl_->siteGroups();
}

bool AtomStereopermutator::thermalized() const {
  return pImpl_->thermalized();
}

const Stereopermutators::Abstract& AtomStereopermutator::getAbstract() const {
  return pImpl_->getAbstract();
}

const Stereopermutators::Feasible& AtomStereopermutator::getFeasible() const {
  return pImpl_->getFeasible();
}

const RankingInformation& AtomStereopermutator::getRanking() const {
  return pImpl_->getRanking();
}

Shapes::Shape AtomStereopermutator::getShape() const {
  return pImpl_->getShape();
}

const AtomStereopermutator::ShapeMap& AtomStereopermutator::getShapePositionMap() const {
  return pImpl_->getShapePositionMap();
}

unsigned AtomStereopermutator::numAssignments() const {
  return pImpl_->numAssignments();
}

unsigned AtomStereopermutator::numStereopermutations() const {
  return pImpl_->numStereopermutations();
}

/* Operators */
bool AtomStereopermutator::operator == (const AtomStereopermutator& other) const {
  return *pImpl_ == *other.pImpl_;
}

bool AtomStereopermutator::operator != (const AtomStereopermutator& other) const {
  return !(*pImpl_ == *other.pImpl_);
}
bool AtomStereopermutator::operator < (const AtomStereopermutator& other) const {
  return *pImpl_ < *other.pImpl_;
}

} // namespace Molassembler
} // namespace Scine
