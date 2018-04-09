#ifndef INCLUDE_TEMPLATE_MAGIC_ADAPTORS_H
#define INCLUDE_TEMPLATE_MAGIC_ADAPTORS_H

#include "Traits.h"
#include "Invoke.h"

/*! @file
 * Boost range like adaptors. Not nearly as good.
 */

namespace temple {

namespace detail {

template<typename Container>
struct SequentialPairGenerator {
  const Container* containerPtr;

  SequentialPairGenerator(const Container& container) : containerPtr(&container) {}

  using ContainerValueType = traits::getValueType<Container>;
  using PairType = std::pair<
    const ContainerValueType&,
    const ContainerValueType&
  >;

  using ContainerIteratorType = decltype(std::begin(std::declval<const Container>()));

  using IteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    PairType, // value_type
    int, // difference_type
    const PairType*, // pointer
    const PairType& // reference
  >;

  class Iterator : public IteratorType {
  private:
    ContainerIteratorType _left, _right;

  public:
    Iterator() = default;
    Iterator(ContainerIteratorType&& left, ContainerIteratorType&& right)
      : _left {left},
        _right {right}
    {}

    void increment() {
      ++_left;
      ++_right;
    }

    Iterator& operator ++ () {
      ++_left;
      ++_right;
      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    /*
    std::enable_if_t<
      std::is_same<
        typename std::iterator_traits<ContainerIteratorType>::iterator_category,
        std::bidirectional_iterator_tag
      >::value,
      Iterator&
    > operator -- () {
      --_left;
      --_right;
      return *this;
    }

    std::enable_if_t<
      std::is_same<
        typename std::iterator_traits<ContainerIteratorType>::iterator_category,
        std::bidirectional_iterator_tag
      >::value,
      Iterator
    > operator -- (int) {
      Iterator prior = *this;
      --(*this);
      return prior;
    }*/

    bool operator == (const Iterator& other) const {
      return _left == other._left && _right == other._right;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {
        *_left,
        *_right
      };
    }
  };

  Iterator begin() const {
    auto second = ++containerPtr->begin();

    return {
      containerPtr->begin(),
      std::move(second)
    };
  }

  Iterator end() const {
    auto prior = --containerPtr->end();

    return {
      std::move(prior),
      containerPtr->end()
    };
  }
};

template<typename Container>
struct AllPairsGenerator {
  const Container* containerPtr;

  AllPairsGenerator(const Container& container) : containerPtr(&container) {}

  using ContainerValueType = traits::getValueType<Container>;
  using PairType = std::pair<
    const ContainerValueType&,
    const ContainerValueType&
  >;

  using ContainerIteratorType = decltype(std::begin(std::declval<const Container>()));

  using IteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    PairType, // value_type
    int, // difference_type
    const PairType*, // pointer
    const PairType& // reference
  >;

  template<typename ContainerIterator>
  class Iterator : public IteratorType {
  private:
    ContainerIterator _left, _right, _end;

  public:
    Iterator() = default;
    Iterator(ContainerIterator&& left, ContainerIterator&& right, ContainerIterator&& end)
      : _left {left},
        _right {right},
        _end {end}
    {}

    Iterator& operator ++ () {
      ++_right;
      if(_right == _end) {
        ++_left;
        _right = _left;
        ++_right;
      }

      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return _left == other._left && _right == other._right;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {
        *_left,
        *_right
      };
    }
  };

  Iterator<ContainerIteratorType> begin() const {
    auto second = ++containerPtr->begin();

    return {
      containerPtr->begin(),
      std::move(second),
      containerPtr->end()
    };
  }

  Iterator<ContainerIteratorType> end() const {
    auto prior = --containerPtr->end();

    return {
      std::move(prior),
      containerPtr->end(),
      containerPtr->end()
    };
  }
};

template<typename ContainerT, typename ContainerU>
struct Zipper {
  const ContainerT* containerTPtr;
  const ContainerU* containerUPtr;

  Zipper(
    const ContainerT& containerT,
    const ContainerU& containerU
  ) : containerTPtr(&containerT),
      containerUPtr(&containerU)
  {}

  using T = traits::getValueType<ContainerT>;
  using U = traits::getValueType<ContainerU>;

  using PairType = std::pair<
    const T&,
    const U&
  >;

  using ContainerTIteratorType = decltype(std::begin(std::declval<const ContainerT>()));
  using ContainerUIteratorType = decltype(std::begin(std::declval<const ContainerU>()));

  using IteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    PairType, // value_type
    int, // difference_type
    const PairType*, // pointer
    const PairType& // reference
  >;

  class Iterator : public IteratorType {
  private:
    ContainerTIteratorType _left;
    ContainerUIteratorType _right;

  public:
    Iterator() = default;
    Iterator(ContainerTIteratorType&& left, ContainerUIteratorType&& right)
      : _left {left},
        _right {right}
    {}

    Iterator& operator ++ () {
      ++_left;
      ++_right;

      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return _left == other._left && _right == other._right;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    PairType operator * () const {
      return {
        *_left,
        *_right
      };
    }
  };

  Iterator begin() const {
    return {
      containerTPtr->begin(),
      containerUPtr->begin()
    };
  }

  Iterator end() const {
    auto tSize = std::distance(containerTPtr->begin(), containerTPtr->end());
    auto uSize = std::distance(containerUPtr->begin(), containerUPtr->end());

    if(tSize == uSize) {
      return {
        containerTPtr->end(),
        containerUPtr->end()
      };
    }

    if(tSize < uSize) {
      auto uIter = containerUPtr->begin();
      std::advance(uIter, tSize);

      return {
        containerTPtr->end(),
        std::move(uIter)
      };
    }

    auto tIter = containerTPtr->begin();
    std::advance(tIter, uSize);

    return {
      std::move(tIter),
      containerUPtr->end()
    };
  }
};

template<typename Container, typename UnaryFunction>
struct Transformer {
  const Container* containerPtr;
  UnaryFunction function;

  Transformer(
    const Container& container,
    UnaryFunction&& function
  ) : containerPtr(&container), function(function) {}

  using ContainerValueType = traits::getValueType<Container>;
  using ReturnType = decltype(
    invoke(
      function,
      std::declval<ContainerValueType>()
    )
  );

  using ContainerIteratorType = decltype(std::begin(std::declval<const Container>()));

  using IteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    ReturnType, // value_type
    int, // difference_type
    const ReturnType*, // pointer
    const ReturnType& // reference
  >;

  class Iterator : public IteratorType {
  private:
    const Transformer* _basePtr;
    ContainerIteratorType _iter;

  public:
    Iterator() = default;
    Iterator(
      const Transformer& base,
      ContainerIteratorType&& iter
    ) : _basePtr(&base),
        _iter {iter}
    {}

    Iterator& operator ++ () {
      ++_iter;
      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    bool operator == (const Iterator& other) const {
      return _iter == other._iter;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    ReturnType operator * () const {
      return invoke(
        _basePtr->function,
        *_iter
      );
    }
  };

  Iterator begin() const {
    return {
      *this,
      containerPtr->begin()
    };
  }

  Iterator end() const {
    return {
      *this,
      containerPtr->end()
    };
  }
};


} // namespace detail

template<typename Container>
auto sequentialPairs(const Container& container) {
  return detail::SequentialPairGenerator<Container>(container);
}

template<typename Container>
auto allPairs(const Container& container) {
  return detail::AllPairsGenerator<Container>(container);
}

template<typename ContainerT, typename ContainerU>
auto zip(const ContainerT& containerT, const ContainerU& containerU) {
  return detail::Zipper<ContainerT, ContainerU>(containerT, containerU);
}

template<typename Container, typename UnaryFunction>
auto transform(
  const Container& container,
  UnaryFunction&& function
) {
  return detail::Transformer<Container, UnaryFunction>(
    container,
    std::forward<UnaryFunction>(function)
  );
}

} // namespace temple

#endif
