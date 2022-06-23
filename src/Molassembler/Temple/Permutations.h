/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides functionality related to permutations.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_PERMUTATIONS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_PERMUTATIONS_H

#include "Molassembler/Temple/Traits.h"
#include "Molassembler/Temple/OperatorSuppliers.h"
#include <algorithm>
#include <tuple>
#include <stdexcept>

namespace Scine {
namespace Molassembler {
namespace Temple {

/*! @brief Calculate the index of permutation of elements in a container
 *
 * The elements in the container must be comparable with a weak ordering.
 *
 * @complexity{@math{\Theta(N^2)}}
 * @note Requires that Container implements operator[](U), where U is implicitly
 *   convertible from std::size_t.
 */
template<class Container>
constexpr std::size_t permutationIndex(const Container& container) {
  const std::size_t size = container.size();

  std::size_t index = 0;
  std::size_t position = 2;// position 1 is paired with factor 0 and so is skipped
  std::size_t factor = 1;

  for(std::size_t p = size - 2; p != std::numeric_limits<std::size_t>::max(); --p) {
    std::size_t largerSuccessors = 0;

    for(std::size_t q = p + 1; q < size; ++q) {
      if(container[q] < container[p]) {
        ++largerSuccessors;
      }
    }

    index += (largerSuccessors * factor);
    factor *= position;
    ++position;
  }

  return index;
}

/*! @brief Index-based in-place swapping of elements in an array-like container
 *
 * @complexity{@math{\Theta(1)}}
 */
template<typename Container>
constexpr void inPlaceSwap(
  Container& data,
  const std::size_t a,
  const std::size_t b
) {
  auto intermediate = std::move(data.at(b));
  data.at(b) = std::move(data.at(a));
  data.at(a) = std::move(intermediate);
}

/*! @brief Index-based in-place reversal of elements in an array-like container
 *
 * @complexity{@math{\Theta(N)}}
 */
template<typename Container>
constexpr void inPlaceReverse(
  Container& data,
  const std::size_t indexFrom,
  const std::size_t indexTo
) {
  std::size_t a = indexFrom;
  std::size_t b = indexTo;
  while(a != b && a != --b) {
    inPlaceSwap(data, a++, b);
  }
}

/*! @brief In-place next permutation
 *
 * In-place next permutation of elements in an array-like type. 1:1 index-based
 * variant of std::next_permutation
 * NOTE: works with std::array only in C++17 (missing constexpr markers)
 *
 * @complexity{@math{O(N/2)}}
 */
template<typename Container>
constexpr bool inPlaceNextPermutation(
  Container& data,
  const std::size_t first,
  const std::size_t last
) {
  std::size_t size = data.size();
  if(!(
    first < last
    && first < size
    && last <= size
  )) {
    throw "Call parameters to inPlaceNextPermutation make no sense!";
  }

  std::size_t i = last - 1;
  std::size_t j = 0;
  std::size_t k = 0;

  while(true) {
    j = i;

    if(
      i != 0
      && data.at(--i) < data.at(j)
    ) {
      k = last;

      while(
        k != 0
        && !(data.at(i) < data.at(--k))
      ) {}

      inPlaceSwap(data, i, k);
      inPlaceReverse(data, j, last);
      return true;
    }

    if(i == first) {
      inPlaceReverse(data, first, last);
      return false;
    }
  }
}

//! @overload
template<typename Container>
constexpr bool inPlaceNextPermutation(Container& data) {
  return inPlaceNextPermutation(data, 0, data.size());
}

/*! @brief In-place previous permutation
 *
 * In-place previous permutation of elements in an array-like type. 1:1
 * index-based variant of std::prev_permutation
 * NOTE: works with std::array only in C++17
 *
 * @complexity{@math{O(N/2)}}
 */
template<typename Container>
constexpr bool inPlacePreviousPermutation(
  Container& data,
  const std::size_t first,
  const std::size_t last
) {
  unsigned size = data.size();
  if(!(
    first < last
    && first < size
    && last <= size
  )) {
    throw "Call parameters to inPlaceNextPermutation make no sense!";
  }

  std::size_t i = last - 1;
  std::size_t j = 0;
  std::size_t k = 0;

  while(true) {
    j = i;

    if(
      i != 0
      && data.at(j) < data.at(--i)
    ) {
      k = last;

      while(
        k != 0
        && !(data.at(--k) < data.at(i))
      ) {}

      inPlaceSwap(data, i, k);
      inPlaceReverse(data, j, last);
      return true;
    }

    if(i == first) {
      inPlaceReverse(data, first, last);
      return false;
    }
  }
}

//! @overload
template<typename Container>
constexpr bool inPlacePreviousPermutation(Container& data) {
  return inPlacePreviousPermutation(data, 0, data.size());
}


//! Calls std::next_permutation
template<class Container>
bool next_permutation(Container& container) {
  return std::next_permutation(
    std::begin(container),
    std::end(container)
  );
}

//! Calls std::prev_permutation
template<class Container>
bool prev_permutation(Container& container) {
  return std::prev_permutation(
    std::begin(container),
    std::end(container)
  );
}

/*!
 * @brief For when you have to implement variable-depth for loops, each with
 *   different limits.
 *
 * Increments a container containing indices for each depth according to a const
 * container containing the limits for each. Limits are exclusive!
 *
 * E.g. for the limits 232, increments through the sequence
 * 000 -> 001 -> 010 -> 011 -> 020 -> 021 -> 100 -> 101 -> 110 -> 111 -> 120 ...
 *
 * @note Requires that Container implements operator[](U), where U is implicitly
 *   convertible from unsigned.
 *
 * @complexity{@math{\Theta(N)}}
 */
template<class Container>
bool nextCombinationPermutation(
  Container& toPermute,
  const Container& limits
) {
  assert(toPermute.size() == limits.size());
  const unsigned cols = toPermute.size();

  // Check if all columns are full
  bool allFull = true;
  for(unsigned i = 0; i < cols; ++i) {
    if(toPermute[i] != limits[i]) {
      allFull = false;
      break;
    }
  }

  if(allFull) {
    return false;
  }

  // Make next permutation
  for(int i = cols - 1; i >= 0; --i) {
    if(toPermute[i] == limits[i]) {
      toPermute[i] = 0;
    } else {
      ++toPermute[i];
      return true;
    }
  }

  return true;
}

/**
 * @brief Container-abstracted permutation
 *
 * Permutations contain only non-negative integers, with each number up less
 * than its length used once. This class can be used only to represent
 * permutations of sequences of a set of items, not of a multiset (with
 * repeated elements).
 *
 * @tparam Container implementing begin/end, at and size methods. If these
 * methods of the container are constexpr, then this class is constexpr.
 *
 * The permutation itself represented in one-line 'notation'. E.g. a vector
 * containing the indices 0, 3, 2, 1 represents the permutation with p(0) = 0,
 * p(1) = 3, p(2) = 2 and p(3) = 1.
 */
struct Permutation : public Crtp::LexicographicComparable<Permutation> {
  using Sigma = std::vector<unsigned>;

  Permutation() = default;

  // Construct from data (unchecked)
  explicit Permutation(Sigma p) : sigma(std::move(p)) {}

  template<typename T, std::enable_if_t<std::is_convertible<T, unsigned>::value, int>* = nullptr>
  explicit Permutation(std::initializer_list<T> vs) {
    sigma.reserve(vs.size());
    for(auto v : vs) {
      sigma.push_back(v);
    }
  }

  //! Construct the i-th permutation of size N
  Permutation(const std::size_t N, std::size_t i) : sigma(N) {
    std::vector<std::size_t> factorials(N);
    factorials.at(0) = 1;
    for(std::size_t k = 1; k < N; ++k) {
      factorials.at(k) = factorials.at(k - 1) * k;
    }

    for(std::size_t k = 0; k < N; ++k) {
      const std::size_t fac = factorials.at(N - 1 - k);
      sigma.at(k) = i / fac;
      i %= fac;
    }

    for(int k = static_cast<int>(N) - 1; k > 0; --k) {
      for(int j = k - 1; j >= 0; --j) {
        if(sigma.at(j) <= sigma.at(k)) {
          ++sigma.at(k);
        }
      }
    }
  }

  //! Construct the identity permutation of size N
  static Permutation identity(unsigned N) {
    Sigma sigma;
    sigma.resize(N);
    for(unsigned i = 0; i < N; ++i) {
      sigma.at(i) = i;
    }
    return Permutation {std::move(sigma)};
  }

  //! Construct a permutation from data (unchecked!)
  template<typename Container>
  static std::enable_if_t<!std::is_same<Container, Sigma>::value, Permutation> from(const Container& p) {
    static_assert(
      std::is_convertible<Traits::getValueType<Container>, unsigned>::value,
      "Value type of container must be convertible to unsigned!"
    );
    Sigma sigma;
    for(auto v : p) {
      sigma.emplace_back(static_cast<unsigned>(v));
    }
    return Permutation {std::move(sigma)};
  }

  template<typename Container>
  static std::enable_if_t<std::is_same<std::decay_t<Container>, Sigma>::value, Permutation> from(Container&& p) {
    return Permutation {std::forward<Container>(p)};
  }

  //! Discover on ordering permutation of an unordered container
  template<typename OtherContainer>
  static Permutation ordering(const OtherContainer& unordered) {
    Permutation p = Permutation::identity(unordered.size());
    std::sort(
      std::begin(p.sigma),
      std::end(p.sigma),
      [&](const auto i, const auto j) {
        return unordered.at(i) < unordered.at(j);
      }
    );
    return p.inverse();
  }

  bool next() {
    return inPlaceNextPermutation(sigma);
  }

  bool prev() {
    return inPlacePreviousPermutation(sigma);
  }

  unsigned index() const {
    return permutationIndex(sigma);
  }

  unsigned indexOf(unsigned v) const {
    const auto begin = std::begin(sigma);
    const auto end = std::end(sigma);
    const auto findIter = std::find(begin, end, v);

    if(findIter == end) {
      throw std::out_of_range("Item not found");
    }

    return findIter - begin;
  }

  Permutation inverse() const {
    const unsigned size = sigma.size();
    Sigma inv;
    inv.resize(size);
    for(unsigned i = 0; i < size; ++i) {
      inv.at(sigma.at(i)) = i;
    }
    return Permutation {inv};
  }

  auto at(const unsigned i) const -> decltype(auto) {
    return sigma.at(i);
  }

  auto operator() (const unsigned i) const -> decltype(auto) {
    return sigma.at(i);
  }

  void push() {
    sigma.push_back(sigma.size());
  }

  void clear() {
    sigma.clear();
  }

  unsigned size() const {
    return sigma.size();
  }

  /**
   * @brief Applies a permutation to the elements of a container by index
   *
   * @param other The container to be permuted
   *
   * @example
   */
  template<typename OtherContainer>
  OtherContainer apply(const OtherContainer& other) const {
    const unsigned otherSize = other.size();
    if(otherSize > size()) {
      throw std::invalid_argument("Container of elements to apply permutation to larger than permutation itself");
    }
    OtherContainer result(otherSize);
    for(unsigned i = 0; i < otherSize; ++i) {
      result.at(sigma.at(i)) = other.at(i);
    }
    return result;
  }

  /*! @brief Compose two permutations
   *
   * The resulting permutation of this composition applies @p this first, then
   * @p other. Note that composition is not commutative!
   */
  Permutation compose(const Permutation& other) const {
    if(size() != other.size()) {
      throw std::invalid_argument("Mismatched permutation sizes!");
    }
    return Permutation::from(inverse().apply(other.sigma));
  }

  auto tie() const {
    return std::tie(sigma);
  }

  //! The permutation in 'one-line' representation
  Sigma sigma;
};

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
