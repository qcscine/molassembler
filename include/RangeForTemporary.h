#ifndef INCLUDE_RANGE_BASED_FOR_LOOP_TEMPORARY_OBJECT_H
#define INCLUDE_RANGE_BASED_FOR_LOOP_TEMPORARY_OBJECT_H

#include <tuple>

template<typename Iterator>
struct RangeForTemporary {
  Iterator _begin, _end;

  RangeForTemporary(Iterator begin, Iterator end) 
  : _begin(std::move(begin)),
    _end(std::move(end))
  {}

  explicit RangeForTemporary(
    std::pair<Iterator, Iterator> iterators
  ) : _begin(std::move(iterators.first)),
      _end(std::move(iterators.second))
  {}

  Iterator& begin() {
    return _begin;
  }

  Iterator& end() {
    return _end;
  }
};

#endif
