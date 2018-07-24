#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_OPTIONALS_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_OPTIONALS_H

#include <boost/optional.hpp>
#include "temple/Traits.h"

/*! @file
 *
 * Implements optional-returning function composition syntactic sugar to avoid
 * repetitive patterns when dealing with lots of optionals.
 */

namespace temple {

namespace detail {

template<typename T>
using Optional = boost::optional<T>;

} // namespace detail

/* Optional-returning function composition */

struct InjectPlaceholder {};

/*!
 * If you want a position in a function call to be replaced with the value of
 * the previous optional, supply this instance in the function call at that
 * position. E.g.
 * \code{.cpp}
 *
 *   optional<double> safe_sqrt(double x) {
 *     return (x >= 0.0) ? std::sqrt(x) : boost::none;
 *   }
 *
 *   optional<double> safe_reciprocal(double x) {
 *     return (x != 0.0) ? 1 / x : boost::none;
 *   }
 *
 *   safe_sqrt( 4.0) | callIfSome(safe_reciprocal, ANS) // yields Some 0.5
 *   safe_sqrt(-4.0) | callIfSome(safe_reciprocal, ANS) // yields None
 *   safe_sqrt( 0.0) | callIfSome(safe_reciprocal, ANS) // yields None
 *
 * \endcode
 */
constexpr InjectPlaceholder ANS;

/*!
 * Class enabling calling a suspended function if the optional left of a binary
 * operator has a value. Permits propagation of results by passing
 * ANS as the suspended function argument at the position where
 * the previous optional value should be inserted.
 *
 * Do not instantiate this directly, prefer using callIfSome function that
 * deduces the CallIfSome class template parameters from the function arguments.
 */
template<typename PartialFunction, typename ... Parameters>
struct CallIfSome {
  using TupleType = std::tuple<Parameters...>;

  // Function pointer function returns what type?
  using FPtrReturnType = typename traits::FunctionPointerReturnType<PartialFunction>::type;
  // Get value type contained in the return optional
  using ReturnType = typename FPtrReturnType::value_type;

  PartialFunction partialFunction;
  TupleType parameters;

  CallIfSome(
    PartialFunction passPartial,
    Parameters&&... passParameters
  ) : partialFunction(passPartial),
      parameters(passParameters...) {}

  /* This class cannot anticipate what kind of result is stored in the left
   * optional, this has to be deduced when the operator | calls value, so the
   * type R is not necessarily the same as ReturnType!
   */
  template<typename R>
  detail::Optional<ReturnType> value(const R& previousResult) const {
    return injectedCallHelper(
      previousResult,
      std::make_index_sequence<
        std::tuple_size<TupleType>::value
      >()
    );
  }

  template<typename R, typename T>
  std::enable_if_t<
    std::is_same<std::decay_t<T>, InjectPlaceholder>::value,
    R
  > replaceIfPlaceholder(
    const R& previousResult,
    const T& /* parameter */
  ) const {
    return previousResult;
  }

  template<typename R, typename T>
  std::enable_if_t<
    !std::is_same<std::decay_t<T>, InjectPlaceholder>::value,
    const T&
  > replaceIfPlaceholder(
    const R& /* previousResult */,
    const T& parameter
  ) const {
    return parameter;
  }

  template<typename R, size_t ... Inds> detail::Optional<ReturnType> injectedCallHelper(
    const R& previousResult,
    std::index_sequence<Inds...>
  ) const {
    return partialFunction(
      replaceIfPlaceholder(
        previousResult,
        std::get<Inds>(parameters)
      )...
    );
  }
};

template<typename PartialFunction, typename ... Parameters>
CallIfSome<PartialFunction, Parameters...> callIfSome (
  PartialFunction partialFunction,
  Parameters&&... parameters
) {
  return {
    partialFunction,
    std::forward<Parameters>(parameters)...
  };
}

/*!
 * Class enabling calling a suspended function if the optional left of a binary
 * operator does not have a value. Propagates an existing Some value if present,
 * foregoing the suspended function call.
 *
 * Do not instantiate this directly, prefer using callIfNone function that
 * deduces the CallIfNone class template parameters from the function arguments.
 */
template<typename PartialFunction, typename ... Parameters>
struct CallIfNone {
  using TupleType = std::tuple<Parameters...>;
  using FPtrReturnType = typename traits::FunctionPointerReturnType<PartialFunction>::type;
  using ReturnType = typename FPtrReturnType::value_type;

  PartialFunction partialFunction;
  TupleType parameters;

  CallIfNone(
    PartialFunction passPartial,
    Parameters&&... passParameters
  ) : partialFunction(passPartial),
      parameters(passParameters...) {}


  detail::Optional<ReturnType> value() const {
    return callHelper(
      std::make_index_sequence<
        std::tuple_size<TupleType>::value
      >()
    );
  }

  template<size_t ... Inds> detail::Optional<ReturnType> callHelper(
    std::index_sequence<Inds...>
  ) const {
    return partialFunction(
      std::get<Inds>(parameters)...
    );
  }
};

template<typename PartialFunction, typename ... Parameters>
CallIfNone<PartialFunction, Parameters...> callIfNone (
  PartialFunction partialFunction,
  Parameters&&... parameters
) {
  return {
    partialFunction,
    std::forward<Parameters>(parameters)...
  };
}

} // namespace temple

/* Global scope operators */
template<typename T, typename PartialFunction, typename ... Parameters>
auto operator | (
  const temple::detail::Optional<T>& valueOptional,
  const temple::CallIfSome<PartialFunction, Parameters...>& other
) {
  using ReturnType = decltype(other.value(valueOptional.value()));

  if(valueOptional) {
    return ReturnType {other.value(valueOptional.value())};
  }

  return ReturnType {boost::none};
}

template<typename T, typename PartialFunction, typename ... Parameters>
auto operator | (
  const temple::detail::Optional<T>& valueOptional,
  const temple::CallIfNone<PartialFunction, Parameters...>& other
) {
  using ReturnType = decltype(other.value());
  using OptionalValueType = typename ReturnType::value_type;
  static_assert(
    std::is_same<
      std::decay_t<T>,
      std::decay_t<OptionalValueType>
    >::value,
    "Your use of callIfNone cannot keep a preceding optional Some value since "
    "the return type of the right function does not match the left optional. "
    "The idea of callIfNone is to try to supply a value if the left is a None, "
    "but to pass on a Some without evaluating the function."
  );

  if(valueOptional) {
    return ReturnType {valueOptional.value()};
  }

  return ReturnType {other.value()};
}

#endif
