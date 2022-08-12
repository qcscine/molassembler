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

AtomStereopermutator::AtomStereopermutator(
    AtomIndex centerAtom,
    Shapes::Shape shape,
    RankingInformation ranking,
    const FeasiblesGenerator& feasibility,
    const ThermalizationPredicate& thermalization,
  const std::vector<std::vector<SiteIndex>>& siteGroups
) : pImpl_(
  std::make_unique<Impl>(
    centerAtom,
    shape,
    std::move(ranking),
    feasibility,
    thermalization,
    siteGroups
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

bool AtomStereopermutator::thermalized(
  AtomIndex centerAtom,
  Shapes::Shape shape,
  const RankingInformation& ranking,
  const Graph& graph
) {
  return Impl::thermalized(centerAtom, shape, ranking, graph);
}

void AtomStereopermutator::assign(boost::optional<unsigned> assignment,
                                  const std::vector<std::vector<SiteIndex>>& siteGroups) {
  pImpl_->assign(std::move(assignment), std::move(siteGroups));
}

void AtomStereopermutator::assignRandom(Random::Engine& engine) {
  pImpl_->assignRandom(engine);
}

void AtomStereopermutator::applyPermutation(const std::vector<AtomIndex>& permutation) {
  pImpl_->applyPermutation(permutation);
}

boost::optional<AtomStereopermutator::ShapeMap> AtomStereopermutator::fit(
  const SiteCentroids& centroids,
  const FeasiblesGenerator& feasibility,
  const ThermalizationPredicate& thermalization
) {
  return pImpl_->fit(centroids, feasibility, thermalization);
}

boost::optional<AtomStereopermutator::ShapeMap> AtomStereopermutator::fit(
  const AngstromPositions& wrapper,
  const FeasiblesGenerator& feasibility,
  const ThermalizationPredicate& thermalization
) {
  return pImpl_->fit(sitePositions(wrapper), feasibility, thermalization);
}

boost::optional<AtomStereopermutator::ShapeMap> AtomStereopermutator::fit(
  const Graph& graph,
  const AngstromPositions& angstromWrapper
) {
  return fit(
    angstromWrapper,
    Stereopermutators::Feasible::Functor(graph),
    thermalizationFunctor(graph)
  );
}

boost::optional<AtomStereopermutator::PropagatedState> AtomStereopermutator::propagate(
  RankingInformation newRanking,
  boost::optional<Shapes::Shape> shapeOption,
  const FeasiblesGenerator& feasibility,
  const ThermalizationPredicate& thermalization
) {
  return pImpl_->propagate(
    std::move(newRanking),
    std::move(shapeOption),
    feasibility,
    thermalization
  );
}

void AtomStereopermutator::propagateVertexRemoval(const AtomIndex removedIndex) {
  pImpl_->propagateVertexRemoval(removedIndex);
}

void AtomStereopermutator::setShape(
  const Shapes::Shape shape,
  const Graph& graph,
  const std::vector<std::vector<SiteIndex>>& siteGroups
) {
  pImpl_->setShape(
    shape,
    Stereopermutators::Feasible::Functor(graph),
    thermalizationFunctor(graph),
    siteGroups
  );
}

void AtomStereopermutator::setShape(
  const Shapes::Shape shape,
  const FeasiblesGenerator& feasibility,
  const ThermalizationPredicate& thermalization,
  const std::vector<std::vector<SiteIndex>>& siteGroups
) {
  pImpl_->setShape(shape, feasibility, thermalization, siteGroups);
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

AtomStereopermutator::SiteCentroids AtomStereopermutator::sitePositions(
  const AngstromPositions& wrapper,
  const std::vector<std::pair<AtomIndex, AtomIndex>>& substitutions
) const {
  return pImpl_->sitePositions(wrapper, substitutions);
}

bool AtomStereopermutator::thermalized() const {
  return pImpl_->thermalized();
}

const Stereopermutators::Abstract& AtomStereopermutator::getAbstract() const {
  return pImpl_->getAbstract();
}

const std::vector<unsigned>& AtomStereopermutator::getFeasible() const {
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
