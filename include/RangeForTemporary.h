#ifndef INCLUDE_RANGE_BASED_FOR_LOOP_TEMPORARY_OBJECT_H
#define INCLUDE_RANGE_BASED_FOR_LOOP_TEMPORARY_OBJECT_H

#include <tuple>

/*! @file
 *
 * Very minor helper class that creates a .begin/.end class for a pair of
 * iterators so that they can be used in a range-for expression.
 */

template<typename Iterator>
struct RangeForTemporary {
  using IteratorPairType = std::pair<Iterator, Iterator>;

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

/*template<typename Iterator>
RangeForTemporary<Iterator> makeRangeForTemporary(
  std::pair<Iterator, Iterator>&& iterators
) {
  return RangeForTemporary<Iterator>(
    std::forward<
      RangeForTemporary<Iterator>::IteratorPairType
    >(iterators)
  );
}*/

#endif
