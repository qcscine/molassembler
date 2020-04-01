/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Stereopermutators/StereopermutatorListImpl.h"

/* The implementation of the const and non-const iterators proceeds in the
 * following steps:
 *
 * - Write the outer implementation of StereopermutatorList::iterator (which is
 *   privately implemented, just like StereopermutatorList)
 * - Since the implementation of the Impl class for the four template arguments
 *   is similar, write a base class that takes care of the common functionality
 * - Implement the iterator Impl classes in terms of the base class
 * - Add instantiations of the iterator class so that the symbols are in the
 *   library when somebody uses their interfaces from the header without a
 *   definition
 *
 * It's disallowed to instantiate the iterator class for any other Permutator
 * types than the ones for which we explicitly instantiate the definitions here,
 * so I don't think the class can be easily misused.
 */

namespace Scine {
namespace molassembler {

/* Iterator facade implementation */
template<typename Permutator>
StereopermutatorList::iterator<Permutator>::iterator(
  iterator&& other
) noexcept = default;

template<typename Permutator>
StereopermutatorList::iterator<Permutator>&
StereopermutatorList::iterator<Permutator>::operator = (
  iterator&& other
) noexcept = default;

template<typename Permutator>
StereopermutatorList::iterator<Permutator>::iterator(
  const iterator& other
) : impl_(std::make_unique<Impl>(*other.impl_)) {}

template<typename Permutator>
StereopermutatorList::iterator<Permutator>&
StereopermutatorList::iterator<Permutator>::operator = (
  const iterator& other
) {
  *impl_ = *other.impl_;
  return *this;
}

template<typename Permutator>
StereopermutatorList::iterator<Permutator>::~iterator() = default;

template<typename Permutator>
StereopermutatorList::iterator<Permutator>::iterator()
  : impl_(std::make_unique<Impl>()) {}

template<typename Permutator>
StereopermutatorList::iterator<Permutator>::iterator(
  ImplType impl,
  bool begin
) : impl_(std::make_unique<Impl>(impl, begin)) {}

template<typename Permutator>
StereopermutatorList::iterator<Permutator>&
StereopermutatorList::iterator<Permutator>::operator ++ () {
  ++(*impl_);
  return *this;
}

template<typename Permutator>
StereopermutatorList::iterator<Permutator>
StereopermutatorList::iterator<Permutator>::operator ++ (int) {
  auto copy = *this;
  ++(*impl_);
  return copy;
}

template<typename Permutator>
typename StereopermutatorList::iterator<Permutator>::reference
StereopermutatorList::iterator<Permutator>::operator * () const {
  return *(*impl_);
}

template<typename Permutator>
bool StereopermutatorList::iterator<Permutator>::operator == (
  const iterator<Permutator>& other
) const {
  return *impl_ == *other.impl_;
}

template<typename Permutator>
bool StereopermutatorList::iterator<Permutator>::operator != (
  const iterator<Permutator>& other
) const {
  return !(*impl_ == *other.impl_);
}

/* Impl base class */
template<typename Iterator>
struct IteratorWrapper {
  Iterator iterator;

  IteratorWrapper() = default;
  explicit IteratorWrapper(Iterator pass) : iterator(std::move(pass)) {}
  virtual ~IteratorWrapper() = default;

  IteratorWrapper& operator ++ () {
    ++iterator;
    return *this;
  }

  bool operator == (const IteratorWrapper& other) const {
    return iterator == other.iterator;
  }

  /* IteratorWrapper operator ++ (int) and operator != are not needed since
   * pImpled iterator wrappers do not call it
   */
};

using AtomMapType = std::unordered_map<AtomIndex, AtomStereopermutator>;
using BondMapType = std::unordered_map<BondIndex, BondStereopermutator, boost::hash<BondIndex>>;

template<>
struct StereopermutatorList::iterator<AtomStereopermutator>::Impl
  : public IteratorWrapper<typename AtomMapType::iterator>
{
  Impl() = default;
  Impl(StereopermutatorList::Impl& impl, bool begin) {
    if(begin) {
      iterator = std::begin(impl.atomStereopermutators);
    } else {
      iterator = std::end(impl.atomStereopermutators);
    }
  }

  AtomStereopermutator& operator * () const {
    return iterator->second;
  }
};

template<>
struct StereopermutatorList::iterator<const AtomStereopermutator>::Impl
  : public IteratorWrapper<typename AtomMapType::const_iterator>
{
  Impl() = default;
  Impl(const StereopermutatorList::Impl& impl, bool begin) {
    if(begin) {
      iterator = std::begin(impl.atomStereopermutators);
    } else {
      iterator = std::end(impl.atomStereopermutators);
    }
  }

  const AtomStereopermutator& operator * () const {
    return iterator->second;
  }
};

template<>
struct StereopermutatorList::iterator<BondStereopermutator>::Impl
  : public IteratorWrapper<typename BondMapType::iterator>
{
  Impl() = default;
  Impl(StereopermutatorList::Impl& impl, bool begin) {
    if(begin) {
      iterator = std::begin(impl.bondStereopermutators);
    } else {
      iterator = std::end(impl.bondStereopermutators);
    }
  }

  BondStereopermutator& operator * () const {
    return iterator->second;
  }
};

template<>
struct StereopermutatorList::iterator<const BondStereopermutator>::Impl
  : public IteratorWrapper<typename BondMapType::const_iterator>
{
  Impl() = default;
  Impl(const StereopermutatorList::Impl& impl, bool begin) {
    if(begin) {
      iterator = std::begin(impl.bondStereopermutators);
    } else {
      iterator = std::end(impl.bondStereopermutators);
    }
  }

  const BondStereopermutator& operator * () const {
    return iterator->second;
  }
};

template class StereopermutatorList::iterator<AtomStereopermutator>;
template class StereopermutatorList::iterator<const AtomStereopermutator>;
template class StereopermutatorList::iterator<BondStereopermutator>;
template class StereopermutatorList::iterator<const BondStereopermutator>;

} // namespace molassembler
} // namespace Scine
