/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/Stereopermutators/StereopermutatorListImpl.h"

#include "boost/optional.hpp"

namespace Scine {
namespace Molassembler {

StereopermutatorList::StereopermutatorList() : impl_(std::make_unique<Impl>()) {}
StereopermutatorList::StereopermutatorList(StereopermutatorList&& other) noexcept = default;
StereopermutatorList& StereopermutatorList::operator = (StereopermutatorList&& other) noexcept = default;
StereopermutatorList::StereopermutatorList(const StereopermutatorList& other) : impl_(std::make_unique<Impl>(*other.impl_)) {}
StereopermutatorList& StereopermutatorList::operator = (const StereopermutatorList& other) {
  *impl_ = *other.impl_;
  return *this;
}
StereopermutatorList::~StereopermutatorList() = default;

AtomStereopermutator& StereopermutatorList::add(
  AtomStereopermutator stereopermutator
) {
  return impl_->add(std::move(stereopermutator));
}

BondStereopermutator& StereopermutatorList::add(
  BondStereopermutator stereopermutator
) {
  return impl_->add(std::move(stereopermutator));
}


//! Apply an index mapping to the list of stereopermutators
void StereopermutatorList::applyPermutation(const std::vector<AtomIndex>& permutation) {
  impl_->applyPermutation(permutation);
}

void StereopermutatorList::clear() {
  impl_->clear();
}

void StereopermutatorList::clearBonds() {
  impl_->clearBonds();
}

void StereopermutatorList::propagateVertexRemoval(const AtomIndex removedIndex) {
  impl_->propagateVertexRemoval(removedIndex);
}

void StereopermutatorList::remove(const AtomIndex index) {
  impl_->remove(index);
}

void StereopermutatorList::remove(const BondIndex& edge) {
  impl_->remove(edge);
}

void StereopermutatorList::try_remove(const AtomIndex index) {
  impl_->try_remove(index);
}

void StereopermutatorList::try_remove(const BondIndex& edge) {
  impl_->try_remove(edge);
}

AtomStereopermutator& StereopermutatorList::at(const AtomIndex index) {
  return impl_->at(index);
}

BondStereopermutator& StereopermutatorList::at(const BondIndex& index) {
  return impl_->at(index);
}

boost::optional<AtomStereopermutator&> StereopermutatorList::option(const AtomIndex index) {
  return impl_->option(index);
}

boost::optional<BondStereopermutator&> StereopermutatorList::option(const BondIndex& edge) {
  return impl_->option(edge);
}

/* Information */
bool StereopermutatorList::empty() const {
  return impl_->empty();
}

unsigned StereopermutatorList::A() const {
  return impl_->A();
}

unsigned StereopermutatorList::B() const {
  return impl_->B();
}

unsigned StereopermutatorList::size() const {
  return impl_->size();
}

const AtomStereopermutator& StereopermutatorList::at(const AtomIndex index) const {
  return static_cast<const Impl*>(impl_.get())->at(index);
}

const BondStereopermutator& StereopermutatorList::at(const BondIndex& index) const {
  return static_cast<const Impl*>(impl_.get())->at(index);
}

boost::optional<const AtomStereopermutator&> StereopermutatorList::option(const AtomIndex index) const {
  // Need to propagate const on the impl pointer to get the right function call
  return static_cast<const Impl*>(impl_.get())->option(index);
}

boost::optional<const BondStereopermutator&> StereopermutatorList::option(const BondIndex& edge) const {
  // Need to propagate const on the impl pointer to get the right function call
  return static_cast<const Impl*>(impl_.get())->option(edge);
}

/* Iterators */
IteratorRange<StereopermutatorList::AtomStereopermutatorIterator>
StereopermutatorList::atomStereopermutators() {
  return {
    AtomStereopermutatorIterator {*impl_, true},
    AtomStereopermutatorIterator {*impl_, false}
  };
}

IteratorRange<StereopermutatorList::AtomStereopermutatorConstIterator>
StereopermutatorList::atomStereopermutators() const {
  return {
    AtomStereopermutatorConstIterator {*impl_, true},
    AtomStereopermutatorConstIterator {*impl_, false}
  };
}

IteratorRange<StereopermutatorList::BondStereopermutatorIterator>
StereopermutatorList::bondStereopermutators() {
  return {
    BondStereopermutatorIterator {*impl_, true},
    BondStereopermutatorIterator {*impl_, false}
  };
}

IteratorRange<StereopermutatorList::BondStereopermutatorConstIterator>
StereopermutatorList::bondStereopermutators() const {
  return {
    BondStereopermutatorConstIterator {*impl_, true},
    BondStereopermutatorConstIterator {*impl_, false}
  };
}

bool StereopermutatorList::hasZeroAssignmentStereopermutators() const {
  return impl_->hasZeroAssignmentStereopermutators();
}

bool StereopermutatorList::hasUnassignedStereopermutators() const {
  return impl_->hasUnassignedStereopermutators();
}

bool StereopermutatorList::compare(
  const StereopermutatorList& other,
  const AtomEnvironmentComponents componentBitmask
) const {
  return impl_->compare(*other.impl_, componentBitmask);
}

/* Operators */
bool StereopermutatorList::operator == (const StereopermutatorList& other) const {
  return *impl_ == *other.impl_;
}

bool StereopermutatorList::operator != (const StereopermutatorList& other) const {
  return !(*impl_ == *other.impl_);
}

} // namespace Molassembler
} // namespace Scine
