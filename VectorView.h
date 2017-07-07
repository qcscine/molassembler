#ifndef INCLUDE_VECTOR_VIEW_H
#define INCLUDE_VECTOR_VIEW_H

#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>

#include "Numeric.h"

/*! @file
 *
 * Without changing the underlying vector, this class permits filtering and/or
 * sorting the data according to custom lambdas and then iterating through the
 * resulting data set in a range-for compatible fashion.
 */

namespace TemplateMagic {

template<typename ValueType>
class VectorView {

private:
  const std::vector<ValueType>& _baseVectorRef;
  std::function<bool(const ValueType&, const ValueType&)> _sortingLambda;
  std::vector<
    std::function<bool(const ValueType&)>
  > _filters;

  std::vector<unsigned> _indexSequence;

  void _recalculateSequence() {
    _indexSequence.resize(_baseVectorRef.size());
    std::iota(
      _indexSequence.begin(),
      _indexSequence.end(),
      0
    );

    if(_sortingLambda) {
      std::sort(
        _indexSequence.begin(),
        _indexSequence.end(),
        [this](const unsigned& a, const unsigned& b) {
          return _sortingLambda(
            _baseVectorRef[a],
            _baseVectorRef[b]
          );
        }
      );
    }

    if(!_filters.empty()) {
      for(const auto& filterFunction : _filters) {
        _indexSequence.erase(
          std::remove_if(
            _indexSequence.begin(),
            _indexSequence.end(),
            [this, filterFunction](const unsigned& a) {
              return filterFunction(
                _baseVectorRef[a]
              );
            }
          ),
          _indexSequence.end()
        );
      }
    }
  }

public:
  using BaseIteratorType = std::iterator<
    std::forward_iterator_tag, // iterator_category
    ValueType,                 // value_type
    unsigned,                  // difference_type
    const ValueType*,          // pointer
    ValueType                  // reference
  >;
  template<typename PointerType> class iterator : public BaseIteratorType {

  private:
    unsigned _idx;
    const std::vector<ValueType>& _vector;
    const std::vector<unsigned>& _sorting;

  public:
    explicit iterator(
      unsigned idx,
      const std::vector<ValueType>& vector,
      const std::vector<unsigned>& sorting
    ) : 
      _idx(idx),
      _vector(vector),
      _sorting(sorting)
    {}

    iterator& operator ++ () { _idx += 1; return *this; }
    iterator operator++ (int) { iterator retval = *this; ++(*this); return retval; }
    bool operator == (iterator other) const { return _idx == other._idx; }
    bool operator != (iterator other) const { return _idx != other._idx; }

    typename BaseIteratorType::reference operator * () const { 
      return _vector.at(_sorting.at(_idx)); 
    }
  };

  VectorView(const std::vector<ValueType>& data) : _baseVectorRef(data) {}

  VectorView& sort(
    std::function<bool(const ValueType&, const ValueType&)> sortLambda
  ) {
    _sortingLambda = sortLambda;
    _recalculateSequence();
    return *this;
  }

  VectorView& filter(std::function<bool(const ValueType&)> filterLambda) {
    _filters.push_back(filterLambda);
    _recalculateSequence();
    return *this;
  }

  VectorView& subset(const std::vector<unsigned>& indices) {
    // Ensure the passed vector contains no out-of-bounds indices
    assert(TemplateMagic::max(indices) < _baseVectorRef.size());

    _indexSequence = indices;
    return *this;
  }

  void resetFilters(const bool update = true) {
    _filters.clear();
    if(update) _recalculateSequence();
  }

  ValueType at(const unsigned& a) const {
    return _baseVectorRef.at(
      _indexSequence.at(a)
    );
  }

  ValueType operator [] (const unsigned& a) const {
    return _baseVectorRef[
      _indexSequence[a]
    ];
  }

  unsigned size() const {
    return _indexSequence.size();
  }

  iterator<const ValueType*> begin() const { 
    return iterator<const ValueType*>(0, _baseVectorRef, _indexSequence); 
  }

  iterator<const ValueType*> end() const { 
    return iterator<const ValueType*>(size(), _baseVectorRef, _indexSequence); 
  }

};

template<typename Container, class FilterFunction>
VectorView<
  traits::getValueType<Container>
> filter(
  const Container& container,
  const FilterFunction&& filterFunction
) {
  auto view = VectorView<
    traits::getValueType<Container>
  > {container};
  view.filter(filterFunction);
  return view;
}

template<typename Container, class SortFunction>
VectorView<
  traits::getValueType<Container>
> sort(
  const Container& container,
  const SortFunction&& filterFunction
) {
  auto view = VectorView<
    traits::getValueType<Container>
  > {container};
  view.sort(filterFunction);
  return view;
}

template<typename Container, typename IndexContainer>
VectorView<
  traits::getValueType<Container>
> subset(
  const Container& container,
  const IndexContainer& indices
) {
  auto view = VectorView<
    traits::getValueType<Container>
  > {container};
  view.subset(indices);
  return view;
}

} // namespace TemplateMagic

#endif
