// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#include "molassembler/Stereocenters/BondStereocenterImpl.h"

namespace molassembler {

constexpr double BondStereocenter::chiralityConstraintTolerance;
constexpr double BondStereocenter::assignmentAcceptanceDihedralThreshold;

/* Impl */

/* Static Impl member functions */

/* BondStereocenter implementation */
BondStereocenter::BondStereocenter(BondStereocenter&& other) noexcept = default;
BondStereocenter& BondStereocenter::operator = (BondStereocenter&& other) noexcept = default;

BondStereocenter::BondStereocenter(const BondStereocenter& other) : _pImpl(
  std::make_unique<Impl>(*other._pImpl)
) {}
BondStereocenter& BondStereocenter::operator = (const BondStereocenter& other) {
  *_pImpl = *other._pImpl;
  return *this;
}
BondStereocenter::~BondStereocenter() = default;

BondStereocenter::BondStereocenter(
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB,
  const BondIndex& edge
) {
  _pImpl = std::make_unique<Impl>(
    stereocenterA,
    stereocenterB,
    edge
  );
}

void BondStereocenter::assign(boost::optional<unsigned> assignment) {
  _pImpl -> assign(std::move(assignment));
}

void BondStereocenter::assignRandom() {
  _pImpl -> assignRandom();
}

void BondStereocenter::fit(
  const AngstromWrapper& angstromWrapper,
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB
) {
  _pImpl -> fit(
    angstromWrapper,
    stereocenterA,
    stereocenterB
  );
}

boost::optional<unsigned> BondStereocenter::assigned() const {
  return _pImpl -> assigned();
}

bool BondStereocenter::hasSameCompositeOrientation(const BondStereocenter& other) const {
  return _pImpl -> hasSameCompositeOrientation(*other._pImpl);
}

boost::optional<unsigned> BondStereocenter::indexOfPermutation() const {
  return _pImpl -> indexOfPermutation();
}

unsigned BondStereocenter::numAssignments() const {
  return _pImpl -> numAssignments();
}

unsigned BondStereocenter::numStereopermutations() const {
  return _pImpl -> numStereopermutations();
}

std::vector<DistanceGeometry::ChiralityConstraint> BondStereocenter::chiralityConstraints(
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB
) const {
  return _pImpl -> chiralityConstraints(
    stereocenterA,
    stereocenterB
  );
}

std::string BondStereocenter::info() const {
  return _pImpl -> info();
}

std::string BondStereocenter::rankInfo() const {
  return _pImpl -> rankInfo();
}

BondIndex BondStereocenter::edge() const {
  return _pImpl -> edge();
}

void BondStereocenter::setModelInformation(
  DistanceGeometry::SpatialModel& model,
  const AtomStereocenter& stereocenterA,
  const AtomStereocenter& stereocenterB,
  double looseningMultiplier
) const {
  _pImpl -> setModelInformation(model, stereocenterA, stereocenterB, looseningMultiplier);
}

bool BondStereocenter::operator == (const BondStereocenter& other) const {
  return *_pImpl == *other._pImpl;
}

bool BondStereocenter::operator != (const BondStereocenter& other) const {
  return !(
    *_pImpl == *other._pImpl
  );
}

} // namespace molassembler
