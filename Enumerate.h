#ifndef INCLUDE_CONTAINER_ENUMERATE_H
#define INCLUDE_CONTAINER_ENUMERATE_H

#include <memory>

namespace enumerate_detail {

template<
  typename T,
  template<typename, typename> class Container
>
class EnumerateTemporary {
public:
  using ContainerType = Container<T, std::allocator<T>>;  
  struct EnumerationStruct {
    const unsigned index;
    const T& value;

    EnumerationStruct(
      const unsigned& index,
      const T& value
    ) : index(index), value(value) {}
  };

private:
  const ContainerType& _containerRef;

public:
  EnumerateTemporary(
    const ContainerType& container
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
    typename ContainerType::const_iterator _it;
    unsigned _index;

  public:
    explicit iterator(
      typename ContainerType::const_iterator it,
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

template<
  typename T,
  template<typename, typename> class Container
>
enumerate_detail::EnumerateTemporary<T, Container> enumerate(
  const Container<T, std::allocator<T>>& container
) {
  return enumerate_detail::EnumerateTemporary<T, Container> (
    container
  );
}

#endif
