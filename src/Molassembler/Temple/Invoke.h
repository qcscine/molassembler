/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Uniform callable invoke from arguments or tuple of arguments
 *
 * Provides uniform C++14 function-like object invoke from arguments or tuple
 * of arguments.
 */

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_INVOKE_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_INVOKE_H

#include <boost/tuple/tuple.hpp>

#include <tuple>
#include <functional>

namespace Scine {
namespace Molassembler {

//! @brief Template shorthands, optimizers and constexpr data types
namespace Temple {
namespace Detail {

template<typename TupleType, typename Function, size_t ... Inds>
auto boostTupleInvokeHelper(
  Function&& function,
  const TupleType& tuple,
  std::index_sequence<Inds...> /* indices */
) {
  return function(
    boost::get<Inds>(tuple)...
  );
}

} // namespace Detail

//! Invokes a function with all values in a given boost tuple
template<
  typename HT,
  typename TT,
  typename Function
> auto invoke(
  Function&& function,
  const boost::tuples::cons<HT, TT>& tuple
) {
  return Detail::boostTupleInvokeHelper(
    function,
    tuple,
    std::make_index_sequence<
      //std::tuple_size<TupleType>::value
      boost::tuples::length<
        boost::tuples::cons<HT, TT>
      >::value
    >()
  );
}

namespace Detail {

template<class> struct sfinae_true : std::true_type {};

/*!
 * @brief Performs invocation of a function by unpacking a tuple-like argument
 *   as arguments and returns the result.
 */
template<typename Function, typename TupleType, size_t ... Inds>
auto invokeHelper(
  Function&& function,
  const TupleType& tuple,
  std::index_sequence<Inds...> /* indices */
) {
  return std::invoke(function, std::get<Inds>(tuple)...);
}

/*!
 * @brief Purely type-based deduction helper for whether a function is callable
 *   by passing an unpacked tuple as arguments. No function body.
 */
struct InvokeTester {
  template<typename Function, typename TupleType, std::size_t ... Inds>
  constexpr auto operator() (
    Function&& function,
    const TupleType& tuple,
    std::index_sequence<Inds...> /* inds */
  ) -> decltype(
    function(
      std::declval<
        std::tuple_element_t<Inds, TupleType>
      >()...
    )
  );
};

/*!
 * @brief Figure out if Function is callable by unpacking TupleType by
 *   forwarding the types and an integer sequence of the tuple's size
 *
 * @note We have to use the InvokeTester instead of deducing the type from
 *   calling invokeHelper because instantiating it and reaching incompatible
 *   types leads to a compilation error instead of a substitution failure.
 */
template<typename Function, typename TupleType>
static auto isTupleCallableTest(int) -> sfinae_true<
  decltype(
    std::declval<InvokeTester>()(
      std::declval<Function>(),
      std::declval<TupleType>(),
      std::make_index_sequence<std::tuple_size<TupleType>::value>()
    )
  )
>;

template<typename Function, typename TupleType, typename... Args>
static auto isTupleCallableTest(long) -> std::false_type;

/*!
 * @brief Returns std::integral_constant<bool> for if Function is callable
 *   by unpacking a single tuple argument supplied in the argument list.
 *
 * @note The integer-long substitution failure trick is explained in Tricks.rst
 */
template<typename Function, typename TupleType, typename... Args>
struct isTupleCallable : decltype(isTupleCallableTest<Function, TupleType, Args...>(0)) {};

} // namespace Detail

/*!
 * @brief If a callable can be called by unpacking a supplied tuple of
 *   arguments, the tuple is called.
 *
 * @note This function is SFINAE-friendly in that if TupleType isn't tuple-like
 *   or the function does not accept the unpacked arguments, this will only
 *   cause a substitution failure and not a compilation error.
 */
template<
  typename Fn,
  typename TupleType,
  std::enable_if_t<Detail::isTupleCallable<Fn, TupleType>::value, int> = 0
> auto invoke(
  Fn&& function,
  const TupleType& tuple
) {
  return Detail::invokeHelper(
    std::forward<Fn>(function),
    tuple,
    std::make_index_sequence<
      std::tuple_size<TupleType>::value
    >()
  );
}

template<
  typename Fn,
  typename... Args,
  std::enable_if_t<
    sizeof...(Args) != 1 || !Detail::isTupleCallable<Fn, Args...>::value,
    int
  > = 0
>
auto invoke(Fn&& function, Args&& ... args) {
  return std::invoke(function, std::forward<Args>(args)...);
}

namespace Detail {

template<typename Functor>
struct Invoker {
  Functor function;

  Invoker(Functor&& passFunction) : function(passFunction) {}

  template<typename TupleType>
  constexpr auto operator() (const TupleType& tuple) const noexcept(noexcept(invoke(function, tuple))) {
    return invoke(function, tuple);
  }
};

} // namespace Detail

template<typename Functor>
auto make_tuple_callable(Functor&& functor) {
  return Detail::Invoker<Functor>(std::forward<Functor>(functor));
}

} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
