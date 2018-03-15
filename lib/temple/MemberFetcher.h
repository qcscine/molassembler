#ifndef INCLUDE_TEMPLATE_MAGIC_MEMBER_FETCHER_H
#define INCLUDE_TEMPLATE_MAGIC_MEMBER_FETCHER_H

#include "Traits.h"

#include <functional>

/*! @file
 *
 * Without changing the underlying container elements, this proxy container
 * permits referential access to class members via a lambda function. This is
 * an attractive alternative to using Container.h's map() function to get 
 * container class element members since this variant avoids copies.
 */

namespace temple {

template<class Container, class FetcherFunction>
class MemberFetcher {
public:
  using ContainerValueType = traits::getValueType<Container>;

  using FetchReturnType = traits::functionReturnType<
    FetcherFunction,
    ContainerValueType
  >;

  using FetcherFunctionDynamicType = std::function<
    FetchReturnType(const ContainerValueType&)
  >;

private:
  const Container& _containerReference;
  const FetcherFunctionDynamicType _fetcher;

public:
  MemberFetcher(
    const Container& container,
    FetcherFunction&& fetcher
  ) :_containerReference(container), _fetcher(fetcher) {}

  using BaseIteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    FetchReturnType,
    unsigned,                  // difference_type
    const FetchReturnType*,          // pointer
    FetchReturnType                  // reference
  >;

  class iterator : public BaseIteratorType {
  private:
    using UnderlyingIteratorType = decltype(
      std::declval<const Container>().begin()
    );

    const FetcherFunctionDynamicType& _fetcher;
    UnderlyingIteratorType _iter;

  public:
    iterator(
      const FetcherFunctionDynamicType& passFetcher,
      UnderlyingIteratorType&& passIter
    ) : _fetcher(passFetcher),
        _iter(passIter)
    {}

    iterator& operator = (const iterator& other) {
      // TODO ensure iterators refer to same base data (unsure if below works)
      /* assert(
        std::addressof(_fetcher.template target<>())
        == std::addressof(other._fetcher.template target<>())
      );*/

      // TODO ensure fetcher functions are identical, maybe via std::function
      // .target member?
      _iter = other._iter;
      return *this;
    }

    iterator& operator ++ () { 
      ++_iter;
      return *this; 
    }

    iterator operator++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval; 
    }

    bool operator == (iterator other) const { return _iter == other._iter; }
    bool operator != (iterator other) const { return _iter != other._iter; }

    typename BaseIteratorType::reference operator * () const { 
      return _fetcher(*_iter);
    }

  };

  iterator begin() const {
    return iterator(
      _fetcher,
      _containerReference.begin()
    );
  }

  iterator end() const {
    return iterator(
      _fetcher,
      _containerReference.end()
    );
  }
};

template<class Container, class FetcherFunction>
MemberFetcher<Container, FetcherFunction> getMember(
  const Container& container,
  FetcherFunction&& fetcher
) {
  return MemberFetcher<Container, FetcherFunction>(
    container,
    std::forward<FetcherFunction>(fetcher)
  );
}

} // namespace temple

#endif
