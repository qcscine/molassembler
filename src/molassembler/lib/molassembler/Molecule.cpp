/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Molecule/MoleculeImpl.h"
#include "molassembler/RankingInformation.h"

namespace Scine {
namespace molassembler {

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
Molecule::Molecule() noexcept : _pImpl(
  std::make_unique<Impl>()
) {}

Molecule::Molecule(Molecule&& other) noexcept = default;
Molecule& Molecule::operator = (Molecule&& rhs) noexcept = default;

Molecule::Molecule(const Molecule& other) : _pImpl(
  std::make_unique<Impl>(*other._pImpl)
) {}
Molecule& Molecule::operator = (const Molecule& rhs) {
  *_pImpl = *rhs._pImpl;
  return *this;
}

Molecule::~Molecule() = default;

Molecule::Molecule(const Utils::ElementType element) noexcept : _pImpl(
  std::make_unique<Impl>(element)
) {}

Molecule::Molecule(
  const Utils::ElementType a,
  const Utils::ElementType b,
  const BondType bondType
) noexcept : _pImpl(
  std::make_unique<Impl>(a, b, bondType)
) {}

Molecule::Molecule(Graph graph) : _pImpl(
  std::make_unique<Impl>(std::move(graph))
) {}

Molecule::Molecule(
  Graph graph,
  const AngstromPositions& positions,
  const boost::optional<
    std::vector<BondIndex>
  >& bondStereopermutatorCandidatesOptional
) : _pImpl(
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
) : _pImpl(
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
  return _pImpl->addAtom(elementType, adjacentTo, bondType);
}

BondIndex Molecule::addBond(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  return _pImpl->addBond(a, b, bondType);
}

void Molecule::applyPermutation(const std::vector<AtomIndex>& permutation) {
  _pImpl->applyPermutation(permutation);
}

void Molecule::assignStereopermutator(
  const AtomIndex a,
  const boost::optional<unsigned>& assignment
) {
  _pImpl->assignStereopermutator(a, assignment);
}

void Molecule::assignStereopermutator(
  const BondIndex& edge,
  const boost::optional<unsigned>& assignment
) {
  _pImpl->assignStereopermutator(edge, assignment);
}

void Molecule::assignStereopermutatorRandomly(
  const AtomIndex a,
  random::Engine& engine
) {
  _pImpl->assignStereopermutatorRandomly(a, engine);
}

void Molecule::assignStereopermutatorRandomly(
  const BondIndex& e,
  random::Engine& engine
) {
  _pImpl->assignStereopermutatorRandomly(e, engine);
}

std::vector<AtomIndex> Molecule::canonicalize(
  const AtomEnvironmentComponents componentBitmask
) {
  return _pImpl->canonicalize(componentBitmask);
}

void Molecule::removeAtom(const AtomIndex a) {
  _pImpl->removeAtom(a);
}

void Molecule::removeBond(
  const AtomIndex a,
  const AtomIndex b
) {
  _pImpl->removeBond(a, b);
}

void Molecule::removeBond(
  const BondIndex& bond
) {
  _pImpl->removeBond(bond.first, bond.second);
}

bool Molecule::setBondType(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  return _pImpl->setBondType(a, b, bondType);
}

void Molecule::setElementType(
  const AtomIndex a,
  const Utils::ElementType elementType
) {
  _pImpl->setElementType(a, elementType);
}

void Molecule::setShapeAtAtom(
  const AtomIndex a,
  const shapes::Shape shape
) {
  _pImpl->setShapeAtAtom(a, shape);
}


/* Information */
boost::optional<AtomEnvironmentComponents> Molecule::canonicalComponents() const {
  return _pImpl->canonicalComponents();
}

boost::optional<shapes::Shape> Molecule::inferShape(
  const AtomIndex index,
  const RankingInformation& ranking
) const {
  return _pImpl->inferShape(index, ranking);
}

std::string Molecule::dumpGraphviz() const {
  return _pImpl->dumpGraphviz();
}

const Graph& Molecule::graph() const {
  return _pImpl->graph();
}

std::size_t Molecule::hash() const {
  return _pImpl->hash();
}

const StereopermutatorList& Molecule::stereopermutators() const {
  return _pImpl->stereopermutators();
}

StereopermutatorList Molecule::inferStereopermutatorsFromPositions(
  const AngstromPositions& angstromWrapper,
  const boost::optional<
    std::vector<BondIndex>
  >& explicitBondStereopermutatorCandidatesOption
) const {
  return _pImpl->inferStereopermutatorsFromPositions(
    angstromWrapper,
    explicitBondStereopermutatorCandidatesOption
  );
}

bool Molecule::canonicalCompare(
  const Molecule& other,
  const AtomEnvironmentComponents componentBitmask
) const {
  return _pImpl->canonicalCompare(*other._pImpl, componentBitmask);
}

bool Molecule::modularCompare(
  const Molecule& other,
  const AtomEnvironmentComponents componentBitmask
) const {
  return _pImpl->modularCompare(*other._pImpl, componentBitmask);
}

std::string Molecule::str() const {
  return _pImpl->str();
}

RankingInformation Molecule::rankPriority(
  const AtomIndex a,
  const std::vector<AtomIndex>& excludeAdjacent,
  const boost::optional<AngstromPositions>& positionsOption
) const {
  return _pImpl->rankPriority(a, excludeAdjacent, positionsOption);
}

/* Operators */
bool Molecule::operator == (const Molecule& other) const {
  return *_pImpl == *other._pImpl;
}

bool Molecule::operator != (const Molecule& other) const {
  return *_pImpl != *other._pImpl;
}


} // namespace molassembler
} // namespace Scine

std::ostream& operator << (
  std::ostream& os,
  const Scine::molassembler::Molecule& molecule
) {
  os << molecule.str();
  return os;
}
