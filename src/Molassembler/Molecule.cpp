/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Molecule/MoleculeImpl.h"
#include "Molassembler/RankingInformation.h"

namespace Scine {
namespace Molassembler {

Utils::AtomCollection Molecule::applyCanonicalizationMap(
  const std::vector<AtomIndex>& canonicalizationIndexMap,
  const Utils::AtomCollection& atomCollection
) {
  return Impl::applyCanonicalizationMap(
    canonicalizationIndexMap,
    atomCollection
  );
}

/* Molecule interface to Impl call forwards */
Molecule::Molecule() noexcept : pImpl_(
  std::make_unique<Impl>()
) {}

Molecule::Molecule(Molecule&& other) noexcept = default;
Molecule& Molecule::operator = (Molecule&& rhs) noexcept = default;

Molecule::Molecule(const Molecule& other) : pImpl_(
  std::make_unique<Impl>(*other.pImpl_)
) {}
Molecule& Molecule::operator = (const Molecule& rhs) {
  *pImpl_ = *rhs.pImpl_;
  return *this;
}

Molecule::~Molecule() = default;

Molecule::Molecule(const Utils::ElementType element) noexcept : pImpl_(
  std::make_unique<Impl>(element)
) {}

Molecule::Molecule(
  const Utils::ElementType a,
  const Utils::ElementType b,
  const BondType bondType
) noexcept : pImpl_(
  std::make_unique<Impl>(a, b, bondType)
) {}

Molecule::Molecule(Graph graph) : pImpl_(
  std::make_unique<Impl>(std::move(graph))
) {}

Molecule::Molecule(
  Graph graph,
  const AngstromPositions& positions,
  const boost::optional<
    std::vector<BondIndex>
  >& bondStereopermutatorCandidatesOptional
) : pImpl_(
  std::make_unique<Impl>(
    std::move(graph),
    positions,
    bondStereopermutatorCandidatesOptional
  )
) {}

Molecule::Molecule(
  Graph graph,
  StereopermutatorList stereopermutators,
  boost::optional<AtomEnvironmentComponents> canonicalComponentsOption
) : pImpl_(
  std::make_unique<Impl>(
    std::move(graph),
    std::move(stereopermutators),
    std::move(canonicalComponentsOption)
  )
) {}

/* Modifiers */
AtomIndex Molecule::addAtom(
  const Utils::ElementType elementType,
  const AtomIndex adjacentTo,
  const BondType bondType
) {
  return pImpl_->addAtom(elementType, adjacentTo, bondType);
}

BondIndex Molecule::addBond(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  return pImpl_->addBond(a, b, bondType);
}

void Molecule::applyPermutation(const std::vector<AtomIndex>& permutation) {
  pImpl_->applyPermutation(permutation);
}

void Molecule::assignStereopermutator(
  const AtomIndex a,
  const boost::optional<unsigned>& assignment
) {
  pImpl_->assignStereopermutator(a, assignment);
}

void Molecule::assignStereopermutator(
  const BondIndex& edge,
  const boost::optional<unsigned>& assignment
) {
  pImpl_->assignStereopermutator(edge, assignment);
}

void Molecule::assignStereopermutatorRandomly(
  const AtomIndex a,
  Random::Engine& engine
) {
  pImpl_->assignStereopermutatorRandomly(a, engine);
}

void Molecule::assignStereopermutatorRandomly(
  const BondIndex& e,
  Random::Engine& engine
) {
  pImpl_->assignStereopermutatorRandomly(e, engine);
}

std::vector<AtomIndex> Molecule::canonicalize(
  const AtomEnvironmentComponents componentBitmask
) {
  return pImpl_->canonicalize(componentBitmask);
}

void Molecule::removeAtom(const AtomIndex a) {
  pImpl_->removeAtom(a);
}

void Molecule::removeBond(
  const AtomIndex a,
  const AtomIndex b
) {
  pImpl_->removeBond(a, b);
}

void Molecule::removeBond(
  const BondIndex& bond
) {
  pImpl_->removeBond(bond.first, bond.second);
}

bool Molecule::setBondType(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  return pImpl_->setBondType(a, b, bondType);
}

void Molecule::setElementType(
  const AtomIndex a,
  const Utils::ElementType elementType
) {
  pImpl_->setElementType(a, elementType);
}

void Molecule::setShapeAtAtom(
  const AtomIndex a,
  const Shapes::Shape shape
) {
  pImpl_->setShapeAtAtom(a, shape);
}


/* Information */
boost::optional<AtomEnvironmentComponents> Molecule::canonicalComponents() const {
  return pImpl_->canonicalComponents();
}

boost::optional<Shapes::Shape> Molecule::inferShape(
  const AtomIndex index,
  const RankingInformation& ranking
) const {
  return pImpl_->inferShape(index, ranking);
}

std::string Molecule::dumpGraphviz() const {
  return pImpl_->dumpGraphviz();
}

const Graph& Molecule::graph() const {
  return pImpl_->graph();
}

std::size_t Molecule::hash() const {
  return pImpl_->hash();
}

const StereopermutatorList& Molecule::stereopermutators() const {
  return pImpl_->stereopermutators();
}

StereopermutatorList Molecule::inferStereopermutatorsFromPositions(
  const AngstromPositions& angstromWrapper,
  const boost::optional<
    std::vector<BondIndex>
  >& explicitBondStereopermutatorCandidatesOption
) const {
  return pImpl_->inferStereopermutatorsFromPositions(
    angstromWrapper,
    explicitBondStereopermutatorCandidatesOption
  );
}

bool Molecule::canonicalCompare(
  const Molecule& other,
  const AtomEnvironmentComponents componentBitmask
) const {
  return pImpl_->canonicalCompare(*other.pImpl_, componentBitmask);
}

boost::optional<std::vector<AtomIndex>> Molecule::modularIsomorphism(
  const Molecule& other,
  const AtomEnvironmentComponents componentBitmask
) const {
  return pImpl_->modularIsomorphism(*other.pImpl_, componentBitmask);
}

std::string Molecule::str() const {
  return pImpl_->str();
}

RankingInformation Molecule::rankPriority(
  const AtomIndex a,
  const std::vector<AtomIndex>& excludeAdjacent,
  const boost::optional<AngstromPositions>& positionsOption
) const {
  return pImpl_->rankPriority(a, excludeAdjacent, positionsOption);
}

/* Operators */
bool Molecule::operator == (const Molecule& other) const {
  return *pImpl_ == *other.pImpl_;
}

bool Molecule::operator != (const Molecule& other) const {
  return *pImpl_ != *other.pImpl_;
}


} // namespace Molassembler
} // namespace Scine

std::ostream& operator << (
  std::ostream& os,
  const Scine::Molassembler::Molecule& molecule
) {
  os << molecule.str();
  return os;
}
