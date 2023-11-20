/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_COMBINATIONS_ADAPTOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_COMBINATIONS_ADAPTOR_H

#include "Molassembler/Temple/ContainerTraits.h"
#include "Molassembler/Temple/Binding.h"

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Adaptors {
namespace Detail {

template<typename Container>
struct Combinations {
//!@name Types
//!@{
  using ValueType = std::decay_t<
    decltype(*std::begin(std::declval<const Container>()))
  >;

  using ValueList = std::vector<ValueType>;

  using Combination = std::vector<ValueType>;

  struct Iterator {
    using iterator_category = std::forward_iterator_tag;
    using value_type = Combination;
    using difference_type = int;
    using pointer = const Combination*;
    using reference = const Combination&;

    Iterator() = default;
    Iterator(const ValueList& container, unsigned R) : pool(container), r(R) {
      const unsigned n = pool.size();
      if(r <= n) {
        // Set up first combination
        indices.resize(r);
        for(unsigned i = 0; i < r; ++i) {
          indices[i] = i;
        }
      } else {
        // Revert to sentinel state
        indices.clear();
        r = 0;
      }
    }

    static Iterator sentinel(const ValueList& container) {
      return Iterator(container, 0);
    }

    Iterator& operator ++ () {
      const unsigned n = pool.size();
      int i = r - 1;
      for(; i >= 0; --i) {
        if(indices[i] != i + n - r) {
          break;
        }
      }
      if(i == -1) {
        // Revert to sentinel state
        indices.clear();
        r = 0;
      } else {
        // Increment
        indices[i] += 1;
        for(unsigned j = i + 1; j < r; ++j) {
          indices[j] = indices[j - 1] + 1;
        }
      }

      return *this;
    }

    Iterator operator ++ (int) {
      Iterator prior = *this;
      ++(*this);
      return prior;
    }

    Combination operator * () const {
      Combination combination(r);
      for(unsigned i = 0; i < r; ++i) {
        combination.at(i) = pool.at(indices.at(i));
      }
      return combination;
    }

    bool operator == (const Iterator& other) const {
      return indices == other.indices && r == other.r;
    }

    bool operator != (const Iterator& other) const {
      return !(*this == other);
    }

    const ValueList& pool;
    std::vector<unsigned> indices;
    unsigned r;
  };
//!@}

  Combinations(const Container& t, unsigned R)
    : values(std::begin(t), std::end(t)), r(R) {}

//!@name Range
//!@{
  Iterator begin() const {
    return {values, r};
  }

  Iterator end() const {
    return Iterator::sentinel(values);
  }
//!@}

//!@name Data
//!@{
  ValueList values;
  unsigned r;
//!@}
};

} // namespace Detail

template<class Container>
auto combinations(Container&& container, unsigned r = 2) {
  return Detail::Combinations<Container>(std::forward<Container>(container), r);
}

} // namespace Adaptors
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
