#ifndef INCLUDE_CONSTEXPR_MAGIC_TUPLETYPE_H
#define INCLUDE_CONSTEXPR_MAGIC_TUPLETYPE_H

#include "UpperTriangularMatrix.h"

/*! @file
 *
 * Provides type-level computations for types enumerated in a tuple.
 */

namespace constable {

namespace TupleType {

namespace detail {

// Helper type needed for deduced return type template enable-if signatures
enum class enabler_t {};

/* This type can be added as a parameter pack with the desired condition in a
 * type signature. If there is no type in the enable_if, this results in an
 * expression error that leads to the signature being ignored in substitution.
 */
template<bool value>
using EnableIf = typename std::enable_if<value, enabler_t>::type;

//! Value variant handler for functions, returns the function call result.
template<typename T, EnableIf<std::is_function<decltype(T::value)>::value>...>
constexpr auto handleValueVariants() {
  return T::value();
}

//! Value variant handler for data members, returns the member itself.
template<typename T, EnableIf<!std::is_function<decltype(T::value)>::value>...>
constexpr auto handleValueVariants() {
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
> constexpr auto unpackHelper(std::index_sequence<I...>) {
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
> constexpr auto mapHelper(std::index_sequence<Inds...>) {
  return makeArray(
    handleValueVariants<
      TemplateFunction<
        std::tuple_element_t<Inds, TupleType>
      >
    >()...
  );
}

/*!
 * Implementation of a special case of mapping in which all possible pairs of
 * types contained in the supplied tuple are passed to a template function and
 * evaluated.
 */
template<
  typename TupleType,
  template<typename ...> class TemplateFunction,
  size_t ... Inds
> constexpr auto mapAllPairsHelper(
  std::index_sequence<Inds...>
) {
  constexpr size_t N = std::tuple_size<TupleType>::value;

  /* Although you might think this is awful, not immediately consuming the
   * generated index pairs, there is basically no elegant way to generate and
   * consume in two places without incurring another function that needs to be
   * instantiated for every pair of indices separately, which would be worse.
   */
  constexpr auto indexPairs = makeArray(
    UpperTriangularMatrixImpl::index_conversion::toDoubleIndex<N>(Inds)...
  );

  return makeArray(
    handleValueVariants<
      TemplateFunction<
        std::tuple_element_t<indexPairs.at(Inds).first, TupleType>,
        std::tuple_element_t<indexPairs.at(Inds).second, TupleType>
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
  size_t ... Inds
> constexpr unsigned countTypeHelper(
  std::index_sequence<Inds...>
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

  for(unsigned i = 0; i < sizeof...(Inds); ++i) {
    sum += trues.at(i);
  }

  return sum;
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

/*! 
 * Tests whether all types in the tuple return true when evaluated against a
 * predicate
 *
 * TODO this is wrong by construction, allOf should be implemented at array-level
 * in Containers.h, not take a Predicate
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

} // namespace TupleType

} // namespace constable

#endif
