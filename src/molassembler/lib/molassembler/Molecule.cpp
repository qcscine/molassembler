// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Molecule/MoleculeImpl.h"

namespace molassembler {

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

Molecule::Molecule(
  const Delib::ElementType a,
  const Delib::ElementType b,
  const BondType bondType
) noexcept : _pImpl(
  std::make_unique<Impl>(a, b, bondType)
) {}

Molecule::Molecule(OuterGraph graph) : _pImpl(
  std::make_unique<Impl>(std::move(graph))
) {}

Molecule::Molecule(
  OuterGraph graph,
  const AngstromWrapper& positions,
  const boost::optional<
    std::vector<BondIndex>
  >& bondStereocenterCandidatesOptional
) : _pImpl(
  std::make_unique<Impl>(
    std::move(graph),
    positions,
    bondStereocenterCandidatesOptional
  )
) {}

Molecule::Molecule(
  OuterGraph graph,
  StereocenterList stereocenters
) : _pImpl(
  std::make_unique<Impl>(
    std::move(graph),
    std::move(stereocenters)
  )
) {}

/* Modifiers */
AtomIndex Molecule::addAtom(
  const Delib::ElementType elementType,
  const AtomIndex adjacentTo,
  const BondType bondType
) {
  return _pImpl->addAtom(elementType, adjacentTo, bondType);
}

void Molecule::addBond(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  _pImpl->addBond(a, b, bondType);
}

void Molecule::assignStereocenter(
  const AtomIndex a,
  const boost::optional<unsigned>& assignment
) {
  _pImpl->assignStereocenter(a, assignment);
}

void Molecule::assignStereocenter(
  const BondIndex& edge,
  const boost::optional<unsigned>& assignment
) {
  _pImpl->assignStereocenter(edge, assignment);
}

void Molecule::assignStereocenterRandomly(const AtomIndex a) {
  _pImpl->assignStereocenterRandomly(a);
}

void Molecule::assignStereocenterRandomly(const BondIndex& e) {
  _pImpl->assignStereocenterRandomly(e);
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

bool Molecule::setBondType(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  return _pImpl->setBondType(a, b, bondType);
}

void Molecule::setElementType(
  const AtomIndex a,
  const Delib::ElementType elementType
) {
  _pImpl->setElementType(a, elementType);
}

void Molecule::setGeometryAtAtom(
  const AtomIndex a,
  const Symmetry::Name symmetryName
) {
  _pImpl->setGeometryAtAtom(a, symmetryName);
}


/* Information */
Symmetry::Name Molecule::determineLocalGeometry(
  const AtomIndex index,
  const RankingInformation& ranking
) const {
  return _pImpl->determineLocalGeometry(index, ranking);
}

std::string Molecule::dumpGraphviz() const {
  return _pImpl->dumpGraphviz();
}

const OuterGraph& Molecule::graph() const {
  return _pImpl->graph();
}

const StereocenterList& Molecule::stereocenters() const {
  return _pImpl->stereocenters();
}

StereocenterList Molecule::inferStereocentersFromPositions(
  const AngstromWrapper& angstromWrapper,
  const boost::optional<
    std::vector<BondIndex>
  >& explicitBondStereocenterCandidatesOption
) const {
  return _pImpl->inferStereocentersFromPositions(
    angstromWrapper,
    explicitBondStereocenterCandidatesOption
  );
}

bool Molecule::modularCompare(
  const Molecule& other,
  const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
) const {
  return _pImpl->modularCompare(*other._pImpl, comparisonBitmask);
}

RankingInformation Molecule::rankPriority(
  const AtomIndex a,
  const std::vector<AtomIndex>& excludeAdjacent,
  const boost::optional<AngstromWrapper>& positionsOption
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

std::ostream& operator << (
  std::ostream& os,
  const molassembler::Molecule& molecule
) {
  const auto& stereocenters = molecule.stereocenters();

  if(!stereocenters.empty()) {
    os << "Stereocenter information:\n";

    for(const auto& stereocenter : stereocenters.atomStereocenters()) {
      os << stereocenter.info() << "\n";
    }

    for(const auto& stereocenter : stereocenters.bondStereocenters()) {
      os << stereocenter.info() << "\n";
    }
  }

  return os;
}
