#ifndef INCLUDE_TEMPLATE_MAGIC_OPTIONALS_H
#define INCLUDE_TEMPLATE_MAGIC_OPTIONALS_H

#include <boost/optional.hpp>
#include "Traits.h"

namespace TemplateMagic {

namespace detail {

template<typename T>
using Optional = boost::optional<T>;

} // namespace detail

/* Optional-returning function composition */

struct InjectPlaceholder {};

constexpr InjectPlaceholder ANS;

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

  detail::Optional<ReturnType> value(const ReturnType& previousResult) const {
    return injectedCallHelper(
      previousResult,
      std::make_index_sequence<
        std::tuple_size<TupleType>::value
      >()
    );
  }

  template<typename T> 
  std::enable_if_t<
    std::is_same<std::decay_t<T>, InjectPlaceholder>::value,
    ReturnType
  > replaceIfPlaceholder(
    const ReturnType& previousResult,
    const T& parameter __attribute__((unused))
  ) const {
    return previousResult;
  }

  template<typename T> 
  std::enable_if_t<
    !std::is_same<std::decay_t<T>, InjectPlaceholder>::value,
    const T&
  > replaceIfPlaceholder(
    const ReturnType& previousResult __attribute__((unused)),
    const T& parameter
  ) const {
    return parameter;
  }

  template<size_t ... Inds> detail::Optional<ReturnType> injectedCallHelper(
    const ReturnType& previousResult,
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

template<typename PartialFunction, typename ... Parameters>
struct CallIfNone {
  using TupleType = std::tuple<Parameters...>;
  using ReturnType = traits::functionReturnType<PartialFunction, Parameters...>;

  PartialFunction partialFunction;
  TupleType parameters;

  CallIfNone(
    PartialFunction passPartial,
    Parameters&&... passParameters
  ) : partialFunction(passPartial),
      parameters(passParameters...) {}


  ReturnType value() const {
    return callHelper(
      std::make_index_sequence<
        std::tuple_size<TupleType>::value
      >()
    );
  }

  template<size_t ... Inds> detail::Optional<double> callHelper(
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

} // namespace TemplateMagic

// Global scope operator injection
template<typename T, typename PartialFunction, typename ... Parameters>
TemplateMagic::detail::Optional<T> operator | (
  const TemplateMagic::detail::Optional<T>& valueOptional,
  const TemplateMagic::CallIfSome<PartialFunction, Parameters...>& other
) {
  if(valueOptional) {
    return other.value(valueOptional.value());
  }

  return boost::none;
}

template<typename T, typename PartialFunction, typename ... Parameters>
TemplateMagic::detail::Optional<T> operator | (
  const TemplateMagic::detail::Optional<T>& valueOptional,
  const TemplateMagic::CallIfNone<PartialFunction, Parameters...>& other
) {
  if(valueOptional) {
    return valueOptional;
  }

  return other.value();
}

#endif
