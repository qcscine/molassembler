#ifndef INCLUDE_TRAITS_H
#define INCLUDE_TRAITS_H

#include <type_traits>
#include <utility>

namespace Traits {

  /* In lieu of C++17 is_callable, the C++14 solution from:
   * From http://talesofcpp.fusionfenix.com/post-11/true-story-call-me-maybe
   * by Augustín Bergé
   *
   * C++11 and C++03 solutions are also available there if need be.
   */
  namespace detail {
    template <typename T>
    using always_void = void;
   
    template <typename Expr, typename Enable = void>
    struct is_callable_impl
      : std::false_type
    {};
   
    template <typename F, typename ...Args>
    struct is_callable_impl<F(Args...), always_void<std::result_of_t<F(Args...)>>>
      : std::true_type
    {};
  }
   
  template <typename Expr>
  struct is_callable
    : detail::is_callable_impl<Expr>
  {};

}

#endif
