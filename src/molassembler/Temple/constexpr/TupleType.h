/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Provides type-level computations for types enumerated in a tuple.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_TUPLE_TYPE_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_CONSTEXPR_TUPLE_TYPE_H

#include "molassembler/Temple/constexpr/Array.h"

#include <tuple>

namespace Scine {
namespace Temple {

//! @brief template metaprogramming metafunctions on tuple types
namespace Tuples {
namespace detail {

//! Value variant handler for functions, returns the function call result.
template<typename T>
constexpr auto handleValueVariants(
  std::enable_if_t<
    std::is_function<decltype(T::value)>::value,
    int
  >* /* enableIfPtr */ = nullptr
) {
  return T::value();
}

//! Value variant handler for data members, returns the member itself.
template<typename T>
constexpr auto handleValueVariants(
  std::enable_if_t<
    !std::is_function<decltype(T::value)>::value,
    int
  >* /* enableIfPtr */ = nullptr
) {
  return T::value;
}

/*!
 * Implementation of the unpacker that forwards all types in a tuple to a
 * template function that takes all types as parameters. Returns the result
 * of the template function, accepting both value functions and members.
 */
template<
  typename Tuple,
  template<typename ...> class TemplateFunction,
  std::size_t... I
> constexpr auto unpackHelper(std::index_sequence<I...> /* inds */) {
  /* Convert the indices parameter pack into a parameter pack containing all
   * types of the tuple, and return the evaluated template function instantiated
   * with all those types
   */
  return handleValueVariants<
    TemplateFunction<
      std::tuple_element_t<I, Tuple>...
    >
  >();
}

/*!
 * Implementation of the mapper, which evaluates a template function
 * individually for all types contained in the supplied tuple and returns the
 * results collected in a std::array.
 */
template<
  typename TupleType,
  template<typename> class TemplateFunction,
  std::size_t... Inds
> constexpr auto mapHelper(std::index_sequence<Inds...> /* inds */) {
  return makeArray(
    handleValueVariants<
      TemplateFunction<
        std::tuple_element_t<Inds, TupleType>
      >
    >()...
  );
}

/*!
 * Implementation of the counter, which returns how often a type occurs in a
 * supplied tuple type.
 */
template<
  typename TupleType,
  typename T,
  std::size_t ... Inds
> constexpr unsigned countTypeHelper(
  std::index_sequence<Inds...> /* inds */
) {
  constexpr std::array<unsigned, sizeof...(Inds)> trues {{
    static_cast<unsigned>(
      std::is_same<
        std::tuple_element_t<Inds, TupleType>,
        T
      >::value
    )...
  }};

  unsigned sum = 0;

  for(unsigned long i = 0; i < sizeof...(Inds); ++i) {
    sum += trues.at(i);
  }

  return sum;
}

template<typename T>
struct RepeatTypeHelper {
  using BaseType = std::tuple<T>;

  template<std::size_t ... Inds>
  static constexpr auto value(std::index_sequence<Inds...> /* inds */)
  -> std::tuple<std::tuple_element_t<Inds * 0, BaseType>...> {}
};

} // namespace detail

/*!
 * Takes a tuple type and a template function that accepts all of the tuple's
 * contained types at once and unpacks these as template parameters to the
 * template function, returning it's value.
 */
template<
  typename Tuple,
  template<typename ...> class TemplateFunction
> constexpr auto unpackToFunction() {
  return detail::unpackHelper<Tuple, TemplateFunction>(
    std::make_index_sequence<
      std::tuple_size<Tuple>::value
    >()
  );
}

/*!
 * Takes a tuple type and a template function that accepts a single type
 * at a time and returns an array of the return values of the template function
 * called with the tuple types.
 */
template<
  typename TupleType,
  template<typename> class TemplateFunction
> constexpr auto map() {
  return detail::mapHelper<TupleType, TemplateFunction>(
    std::make_index_sequence<
      std::tuple_size<TupleType>::value
    >()
  );
}

//! Counts how often a type is contained in a tuple type.
template<
  typename TupleType,
  typename T
> constexpr unsigned countType() {
  return detail::countTypeHelper<TupleType, T>(
    std::make_index_sequence<
      std::tuple_size<TupleType>::value
    >()
  );
}

/*! @brief all_of with tuple types and template metafunctions
 *
 * Tests whether all types in the tuple return true when evaluated against a
 * predicate
 *
 * @complexity{@math{O(N)}}
 */
template<
  typename TupleType,
  template<typename> class UnaryPredicate
> constexpr bool allOf() {
  constexpr size_t N = std::tuple_size<TupleType>::value;
  constexpr auto mapped = map<TupleType, UnaryPredicate>();
  for(unsigned i = 0; i < N; ++i) {
    if(!mapped.at(i)) {
      return false;
    }
  }

  return true;
}

template<typename T, unsigned repeats>
struct RepeatType {
  using type = decltype(
    detail::RepeatTypeHelper<T>::value(
      std::make_index_sequence<repeats>()
    )
  );
};

} // namespace Tuples
} // namespace Temple
} // namespace Scine

#endif
