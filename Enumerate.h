#ifndef INCLUDE_CONTAINER_ENUMERATE_H
#define INCLUDE_CONTAINER_ENUMERATE_H

#include <memory>

namespace enumerate_detail {

template<class Container>
class EnumerateTemporary {
public:
  // Get the bare Container containing type
  using T = typename std::remove_reference<
    decltype(
      *(
        std::declval<Container>().begin()
      )
    )
  >::type;

  struct EnumerationStruct {
    const unsigned index;
    const T& value;

    EnumerationStruct(
      const unsigned& index,
      const T& value
    ) : index(index), value(value) {}
  };

private:
  const Container& _containerRef;

public:
  EnumerateTemporary(
    const Container& container
  ) : _containerRef(container) 
  {}

  using BaseIteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    EnumerationStruct,         // value_type
    unsigned,                  // difference_type
    const EnumerationStruct*,  // pointer
    const EnumerationStruct    // reference
  >;

  template<typename PointerType> 
  class iterator : public BaseIteratorType {
  private:
    typename Container::const_iterator _it;
    unsigned _index;

  public:
    explicit iterator(
      typename Container::const_iterator it,
      unsigned index
    ) : _it(it), _index(index) {}

    iterator& operator ++ () { 
      ++_it;
      ++_index;
      return *this; 
    }

    iterator operator++ (int) {
      iterator retval = *this;
      ++(*this);
      return retval; 
    }

    bool operator == (iterator other) const {
      return _it == other._it;
    }

    bool operator != (iterator other) const {
      return !(
        *this == other
      );
    }

    typename BaseIteratorType::reference operator * () const { 
      return EnumerationStruct {
        _index,
        *_it
      };
    }
  };

  iterator<const EnumerationStruct*> begin() const {
    return iterator<const EnumerationStruct*>(
      _containerRef.begin(),
      0
    );
  }

  /* consider enable_if-ing this to support containers without .size() by means
   * of an alternate implementation using .end() - .begin()
   */
  iterator<const EnumerationStruct*> end() const {
    return iterator<const EnumerationStruct*>(
      _containerRef.end(),
      _containerRef.size()
    );
  }
};

} // eo namespace enumerate_detail

/*! Returns an EnerateTemporary for use with range-for expressions that
 * generates a struct with members index and value for every contained element.
 * Requires that the container implements begin(), end() and size() members.
 * Should involve minimal copying, mostly uses references, though no space or 
 * time complexity guarantees are given.
 * 
 * We realize this may seem like overkill, particularly when many STL Containers 
 * elements can be accessed with operator [], and customary loops are more
 * adequate. However, perhaps some custom containers do not follow STL
 * conventions or do not implement operator [], and for these, use this.
 */
template<class Container>
enumerate_detail::EnumerateTemporary<Container> enumerate(
  const Container& container
) {
  return enumerate_detail::EnumerateTemporary<Container> (
    container
  );
}

#endif
