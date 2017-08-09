#ifndef INCLUDE_CONSTEXPR_MAGIC_TUPLETYPE_H
#define INCLUDE_CONSTEXPR_MAGIC_TUPLETYPE_H

#include "UpperTriangularMatrix.h"

namespace ConstexprMagic {

namespace TupleType {

namespace detail {

template<
  typename Tuple,
  template<typename ...> class TemplateFunction,
  std::size_t... I
> constexpr auto unpackHelper(std::index_sequence<I...>) {
  /* Convert the indices parameter pack into a parameter pack containing all
   * types of the tuple, and return the evaluated template function instantiated
   * with all those types
   */
  return TemplateFunction<
    std::tuple_element_t<I, Tuple>...
  >::value;
}

template<
  typename TupleType,
  template<typename> class TemplateFunction,
  std::size_t... Inds
> constexpr auto mapHelper(std::index_sequence<Inds...>) {
  using ReturnType = decltype(
    TemplateFunction<
      std::tuple_element_t<0, TupleType>
    >::value
  );

  return std::array<ReturnType, sizeof...(Inds)> {{
    TemplateFunction<
      std::tuple_element_t<Inds, TupleType>
    >::value...
  }};
}

template<
  typename TupleType,
  template<typename ...> class TemplateFunction,
  size_t ... Inds
> constexpr auto mapAllPairsHelper(
  std::index_sequence<Inds...>
) {
  constexpr size_t N = std::tuple_size<TupleType>::value;
  constexpr size_t C = sizeof...(Inds);

  constexpr auto indexPairs = std::array<
    std::pair<size_t, size_t>,
    C
  > {{
    UpperTriangularMatrixImpl::index_conversion::toDoubleIndex<N>(Inds)...
  }};

  constexpr auto results = std::array<unsigned, C> {{
    TemplateFunction<
      std::tuple_element_t<indexPairs.at(Inds).first, TupleType>,
      std::tuple_element_t<indexPairs.at(Inds).second, TupleType>
    >::value...
  }};

  return results;
}

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

/*!
 * Takes a tuple type and a template function that accepts pairs of the tuple's
 * contained types. Returns an array containing the result of the mapping for
 * all distinct possible type pairs.
 *
 * With distinct, we mean that pairs of the same types are not created unless a
 * type is explicitly repeated in the tuple type: Say the tuple contains Apple, 
 * Banana and Orange, then the generated pairs are: Apple & Banana, 
 * Apple & Orange and Banana & Orange.
 */
template<
  typename TupleType,
  template<typename ...> class TemplateFunction
> constexpr auto mapAllPairs() {
  constexpr size_t N = std::tuple_size<TupleType>::value;

  return detail::mapAllPairsHelper<TupleType, TemplateFunction>(
    std::make_index_sequence<(N * N - N) / 2>()
  );
}

} // namespace TupleType

} // namespace ConstexprMagic

#endif
