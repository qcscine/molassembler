/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Strongly typed index permutations
 */
#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_STRONG_INDEX_PERMUTATION_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_STRONG_INDEX_PERMUTATION_H

#include "Molassembler/Temple/Permutations.h"

namespace Scine {
namespace Molassembler {
namespace Temple {

template<typename Key, typename Value>
struct StrongIndexPermutation : public Crtp::LexicographicComparable<StrongIndexPermutation<Key, Value>> {
  //!@name Types
  //!@{
  using key_type = Key;
  using value_type = Value;

  struct Iterator : Crtp::LexicographicComparable<Iterator> {
    using base_iterator = std::vector<unsigned>::const_iterator;
    using iterator_category = base_iterator::iterator_category;
    using value_type = std::pair<Key, Value>;

    Iterator(base_iterator a, base_iterator b) : begin(std::move(a)), iter(std::move(b)) {}

    value_type operator * () const {
      return {
        Key {static_cast<typename Key::value_type>(iter - begin)},
        Value {*iter}
      };
    }

    Iterator& operator ++ () {
      ++iter;
      return *this;
    }

    Iterator operator ++ (int) {
      Iterator copy = *this;
      ++(*this);
      return copy;
    }

    auto tie() const {
      return std::tie(iter);
    }

    base_iterator begin;
    base_iterator iter;
  };

  using iterator = Iterator;
  using const_iterator = Iterator;
  //!@}

  //!@name Constructors
  //!@{
  StrongIndexPermutation() = default;

  explicit StrongIndexPermutation(Permutation p) : permutation(std::move(p)) {}

  explicit StrongIndexPermutation(Permutation::Sigma p) : permutation(std::move(p)) {}

  template<typename T, std::enable_if_t<std::is_convertible<T, unsigned>::value, int>* = nullptr>
  explicit StrongIndexPermutation(std::initializer_list<T> vs) : permutation(std::move(vs)) {}
  //!@}

  //!@name Factories
  //!@{
  template<typename Container>
  static std::enable_if_t<!std::is_same<Container, Permutation::Sigma>::value, StrongIndexPermutation> from(const Container& p) {
    static_assert(
      std::is_convertible<Traits::getValueType<Container>, unsigned>::value,
      "Value type of container must be convertible to unsigned!"
    );
    Permutation::Sigma sigma;
    for(auto v : p) {
      sigma.emplace_back(static_cast<unsigned>(v));
    }
    return StrongIndexPermutation {Permutation {std::move(sigma)}};
  }
  //!@}

  //!@name Information
  //!@{
  Value at(const Key i) const {
    return Value {permutation.at(i)};
  }

  Value operator() (const Key i) const {
    return Value {permutation.at(i)};
  }

  Key indexOf(const Value v) const {
    return Key {permutation.indexOf(v)};
  }

  unsigned size() const {
    return permutation.size();
  }

  StrongIndexPermutation<Value, Key> inverse() const {
    return StrongIndexPermutation<Value, Key> {permutation.inverse()};
  }

  template<typename OtherValue>
  StrongIndexPermutation<Key, OtherValue>
  compose(const StrongIndexPermutation<Value, OtherValue>& other) const {
    return StrongIndexPermutation<Key, OtherValue> {permutation.compose(other.permutation)};
  }

  StrongIndexPermutation compose(const StrongIndexPermutation<Value, Value>& other) const {
    return StrongIndexPermutation {permutation.compose(other.permutation)};
  }

  auto tie() const {
    return std::tie(permutation);
  }
  //!@}

  //!@name Modification
  //!@{
  void clear() {
    permutation.clear();
  }
  //!@}

  //!@name Iterators
  //!@{
  Iterator begin() const {
    return Iterator {
      std::begin(permutation.sigma),
      std::begin(permutation.sigma)
    };
  }

  Iterator end() const {
    return Iterator {
      std::begin(permutation.sigma),
      std::end(permutation.sigma)
    };
  }
  //!@}

  Permutation permutation;
};

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
