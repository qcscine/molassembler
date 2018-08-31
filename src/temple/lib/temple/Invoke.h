#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_INVOKE_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_INVOKE_H

#include <boost/tuple/tuple.hpp>

#include <tuple>
#include <functional>

/*! @file
 *
 * Provides uniform C++14 function-like object invoke from arguments or tuple
 * of arguments.
 */

namespace temple {

namespace detail {

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

} // namespace detail

//! Invokes a function with all values in a given boost tuple
template<
  typename HT,
  typename TT,
  typename Function
> auto invoke(
  Function&& function,
  const boost::tuples::cons<HT, TT>& tuple
) {
  return detail::boostTupleInvokeHelper(
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

namespace detail {

template<class> struct sfinae_true : std::true_type {};

template<class Callable, typename ... Args>
static auto isCallableTest(int) -> sfinae_true<
  typename std::result_of<Callable(Args...)>::type
>;

template<class Callable, typename ... Args>
static auto isCallableTest(long) -> std::false_type;

template<class Callable, typename ... Args>
struct isCallable : decltype(isCallableTest<Callable, Args...>(0)) {};

template<class TupleType>
static auto isTupleLikeTest(int) -> sfinae_true<
  std::tuple_element_t<0, TupleType>
>;

template<class TupleType>
static auto isTupleLikeTest(long) -> std::false_type;

template<class TupleType>
struct isTupleLike : decltype(isTupleLikeTest<TupleType>(0)) {};

template<typename Function, typename TupleType, size_t ... Inds>
auto invokeHelper(
  Function&& function,
  const TupleType& tuple,
  std::index_sequence<Inds...> /* indices */
) {
  return function(
    std::get<Inds>(tuple)...
  );
}

template<
  template<typename...> class TemplateFunction,
  typename Function,
  typename ArgumentTuple,
  size_t ... Inds
> auto unpackHelper(std::index_sequence<Inds...> /* indices */) {
  return TemplateFunction<Function, std::tuple_element_t<Inds, ArgumentTuple>...>::value;
}

template<
  template<typename...> class TemplateFunction,
  typename Function,
  typename ArgumentTuple
> auto unpack() {
  return unpackHelper<TemplateFunction, Function, ArgumentTuple>(
    std::make_index_sequence<std::tuple_size<ArgumentTuple>::value>()
  );
}

template<typename Function, typename TupleType>
static auto isTupleCallableTest(int) -> sfinae_true<
    decltype(
      detail::invokeHelper(
        std::declval<Function>(),
        std::declval<TupleType>(),
        std::make_index_sequence<std::tuple_size<TupleType>::value>()
      )
    )
>;

template<typename Function, typename TupleType, typename... Args>
static auto isTupleCallableTest(long) -> std::false_type;

template<typename Function, typename TupleType, typename... Args>
struct isTupleCallable : decltype(isTupleCallableTest<Function, TupleType, Args...>(0)) {};

} // namespace detail

//! Invokes a function with all values in a given tuple
template<
  typename Function,
  typename TupleType,
  std::enable_if_t<
    detail::isTupleCallable<Function, TupleType>::value,
    int
  > = 0
> auto invoke(
  Function&& function,
  const TupleType& tuple
) {
  return detail::invokeHelper(
    std::forward<Function>(function),
    tuple,
    std::make_index_sequence<
      std::tuple_size<TupleType>::value
    >()
  );
}

// Suggested C++17 std::invoke reference implementation from N4169

template<
  typename Fn,
  typename... Args,
  std::enable_if_t<
    !(sizeof...(Args) == 1 && detail::isTupleCallable<Fn, Args...>::value)
    && std::is_member_pointer<std::decay_t<Fn>>{},
    int
  > = 0
> constexpr decltype(auto) invoke(Fn&& f, Args&&... args)
noexcept(noexcept(std::mem_fn(f)(std::forward<Args>(args)...)))
{
  return std::mem_fn(f)(std::forward<Args>(args)...);
}

template<
  typename Fn,
  typename... Args,
  std::enable_if_t<
    !(sizeof...(Args) == 1 && detail::isTupleCallable<Fn, Args...>::value)
    && !std::is_member_pointer<std::decay_t<Fn>>{},
    int
  > = 0
> constexpr decltype(auto) invoke(Fn&& f, Args&&... args)
    noexcept(noexcept(std::forward<Fn>(f)(std::forward<Args>(args)...)))
{
  return std::forward<Fn>(f)(std::forward<Args>(args)...);
}

namespace detail {

template<typename Functor>
struct Invoker {
  Functor function;

  Invoker(Functor&& passFunction) : function(passFunction) {}

  template<typename TupleType>
  constexpr auto operator() (const TupleType& tuple) const noexcept {
    return invoke(function, tuple);
  }
};

} // namespace detail

template<typename Functor>
auto make_tuple_callable(Functor&& functor) {
  return detail::Invoker<Functor>(std::forward<Functor>(functor));
}


} // namespace temple

#endif
