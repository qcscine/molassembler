#ifndef INCLUDE_TEMPLATE_MAGIC_TRAITS_H
#define INCLUDE_TEMPLATE_MAGIC_TRAITS_H

#include <vector>

namespace TemplateMagic {

namespace traits {

/* Is_callable trait
 *
 * In lieu of C++17 is_callable, the C++14 solution from:
 * From http://talesofcpp.fusionfenix.com/post-11/true-story-call-me-maybe
 * by Augustín Bergé
 *
 * C++11 and C++03 solutions are also available there if need be.
 */
namespace detail {
  template <typename T>
  using always_void = void;
 
  template <typename Expr, typename Enable = void>
  struct is_callable_impl : std::false_type
  {};
 
  template <typename F, typename ...Args>
  struct is_callable_impl<F(Args...), always_void<std::result_of_t<F(Args...)>>>
    : std::true_type
  {};
 
  template <typename Expr>
  struct is_callable : is_callable_impl<Expr>
  {};
} // namespace detail

template<typename Expr>
constexpr bool isCallableValue = detail::is_callable<Expr>::value;

/* Get the base type a container holds via the begin iterator
 */
namespace detail {
  template<class ContainerType> 
  struct getValueTypeImpl {
    using type = typename std::remove_const<
      typename std::remove_reference<
        decltype(
          *(
            std::declval<ContainerType>()
          ).begin()
        )
      >::type
    >::type;
  };

  // Specialization for std::vector<bool>, which returns an awkward proxy object
  template<>
  struct getValueTypeImpl<std::vector<bool>> {
    using type = bool;
  };
} // namespace detail

template<class ContainerType>
using getValueType = typename detail::getValueTypeImpl<ContainerType>::type;

template<class Function, typename ...Args>
using functionReturnType = std::result_of_t<Function(Args...)>;

} // namespace traits

} // namespace TemplateMagic

#endif
