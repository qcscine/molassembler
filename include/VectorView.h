#include <vector>
#include <functional>
#include <algorithm>

template<typename ValueType>
class VectorView {

private:
  std::vector<ValueType>& _baseVectorRef;
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

    if(_filters.size() > 0) {
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
    std::vector<ValueType>& _vector;
    std::vector<unsigned>& _sorting;

  public:
    explicit iterator(
      unsigned idx,
      std::vector<ValueType>& vector,
      std::vector<unsigned>& sorting
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
      return _vector[_sorting[_idx]]; 
    }
  };

  VectorView(std::vector<ValueType>& data) : _baseVectorRef(data) {}

  VectorView sort(
    std::function<bool(const ValueType&, const ValueType&)> sortLambda
  ) {
    _sortingLambda = sortLambda;
    _recalculateSequence();
    return *this;
  }

  VectorView filter(std::function<bool(const ValueType&)> filterLambda) {
    _filters.push_back(filterLambda);
    _recalculateSequence();
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

  iterator<const ValueType*> begin() { 
    return iterator<const ValueType*>(0, _baseVectorRef, _indexSequence); 
  }

  iterator<const ValueType*> end() { 
    return iterator<const ValueType*>(size(), _baseVectorRef, _indexSequence); 
  }

};
