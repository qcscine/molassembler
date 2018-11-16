// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_VECTOR_VIEW_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_VECTOR_VIEW_H

#include "temple/constexpr/Numeric.h"

#include <algorithm>

/*! @file
 *
 * @brief Nonmodifying vector adaptor enabling filtering & sorting
 *
 * Without changing the underlying vector, this class permits filtering and/or
 * sorting the data according to custom lambdas and then iterating through the
 * resulting data set in a range-for compatible fashion.
 */

namespace temple {

template<typename ValueType>
class VectorView {

private:
  //! A reference to the underlying data
  const std::vector<ValueType>& _baseVectorRef;

  //! Only one lambda can be used to sort the data
  std::function<bool(const ValueType&, const ValueType&)> _sortingLambda;

  //! A list of filter lambdas that remove data if the predicate returns true
  std::vector<
    std::function<bool(const ValueType&)>
  > _filters;

  /*!
   * The current resulting index sequence of elements in the underlying data
   * based on sorting AND filtering that is passed to iterators that progress
   * through this list. This is maintained at any change of sorting or filtering
   * lambdas
   */
  std::vector<unsigned> _indexSequence;

  /*!
   * Recalculates the index sequence if a change has occurred in any of the
   * stored lambdas
   */
  void _recalculateSequence() {
    _indexSequence.resize(_baseVectorRef.size());
    std::iota(
      std::begin(_indexSequence),
      std::end(_indexSequence),
      0
    );

    if(_sortingLambda) {
      std::sort(
        std::begin(_indexSequence),
        std::end(_indexSequence),
        [this](const unsigned a, const unsigned b) {
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
            std::begin(_indexSequence),
            std::end(_indexSequence),
            [this, filterFunction](const unsigned a) {
              return filterFunction(
                _baseVectorRef[a]
              );
            }
          ),
          std::end(_indexSequence)
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
    // Maintain the current position in the sorting sequence
    unsigned _idx;
    // Reference to base data
    const std::vector<ValueType>& _vector;
    // Reference to sorting sequence
    const std::vector<unsigned>& _sorting;

  public:
    iterator(
      unsigned idx,
      const std::vector<ValueType>& vector,
      const std::vector<unsigned>& sorting
    ) : _idx(idx),
        _vector(vector),
        _sorting(sorting)
    {}

    iterator& operator = (const iterator& other) {
      // ensure iterators refer to same instance
      assert(
        std::addressof(_vector) == std::addressof(other._vector)
        && std::addressof(_sorting) == std::addressof(other._sorting)
      );

      _idx = other._idx;

      return *this;
    }

    iterator& operator ++ () { _idx += 1; return *this; }
    iterator operator++ (int) { iterator retval = *this; ++(*this); return retval; }
    bool operator == (iterator other) const { return _idx == other._idx; }
    bool operator != (iterator other) const { return _idx != other._idx; }

    typename BaseIteratorType::reference operator * () const {
      return _vector.at(_sorting.at(_idx));
    }
  };

  VectorView(const std::vector<ValueType>& data) : _baseVectorRef(data) {}

  //! Sorts the current data by a passed binary comparator
  VectorView& sort(
    std::function<bool(const ValueType&, const ValueType&)> sortLambda
  ) {
    _sortingLambda = sortLambda;
    _recalculateSequence();
    return *this;
  }

  /*!
   * Filters out data from the underlying vector if the unary predicate
   * returns true
   */
  VectorView& filter(std::function<bool(const ValueType&)> filterLambda) {
    _filters.push_back(filterLambda);
    _recalculateSequence();
    return *this;
  }

  /*!
   * Allows index-based subsetting of the underlying data if previous evaluation
   * has been performed elsewhere
   */
  VectorView& subset(const std::vector<unsigned>& indices) {
    // Ensure the passed vector contains no out-of-bounds indices
    assert(temple::max(indices) < _baseVectorRef.size());

    _indexSequence = indices;
    return *this;
  }

  //! Removes all filters on the data
  void resetFilters(const bool update = true) {
    _filters.clear();
    if(update) _recalculateSequence();
  }

  ValueType at(const unsigned a) const {
    return _baseVectorRef.at(
      _indexSequence.at(a)
    );
  }

  ValueType operator [] (const unsigned a) const {
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

  std::vector<ValueType> getCopy() const {
    return {
      this->begin(),
      this->end()
    };
  }
};

/*!
 * Provides a VectorView proxy object on a container that filters out some
 * elements if the passed unary predicate returns true
 */
template<class Container, class FilterFunction>
VectorView<
  traits::getValueType<Container>
> view_filter(
  const Container& container,
  const FilterFunction&& filterFunction
) {
  auto view = VectorView<
    traits::getValueType<Container>
  > {container};
  view.filter(filterFunction);
  return view;
}

/*!
 * Provides a VectorView proxy object on a container that sorts the underlying
 * data using a provided binary comparator function.
 */
template<class Container, class SortFunction>
VectorView<
  traits::getValueType<Container>
> view_sort(
  const Container& container,
  const SortFunction&& sortingFunction
) {
  auto view = VectorView<
    traits::getValueType<Container>
  > {container};
  view.sort(sortingFunction);
  return view;
}

/*!
 * Provides a VectorView proxy object on a container that subsets the underlying
 * data based on a provided index list.
 */
template<class Container, typename IndexContainer>
VectorView<
  traits::getValueType<Container>
> view_subset(
  const Container& container,
  const IndexContainer& indices
) {
  auto view = VectorView<
    traits::getValueType<Container>
  > {container};
  view.subset(indices);
  return view;
}

} // namespace temple

#endif
