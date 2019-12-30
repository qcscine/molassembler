# Tricks

Throughout the library, a lot of template trickery is used. Some of it is
recurring, and warrants a central explanation.

## Integer-Long-Reinterpretation

Suppose you have a Container type, and you want to figure out whether the
container implements a size() member. In temple, this is implemented as

```c++
template<class>
struct sfinae_true : std::true_type {};

namespace detail {

template<class Container>
static auto testHasSize(int) -> sfinae_true<
  decltype(
    std::declval<Container>().size()
  )
>;

template<class Container>
static auto testHasSize(long) -> std::false_type;

} // namespace detail

template<class Container>
struct hasSize : decltype(detail::testHasSize<Container>(0)) {};
```

`sfinae_true` exists to be conditionally instantiable. If you try to instantiate
it with a `decltype` template argument type expression that results in an
expression error, instantiation fails. If, instead, the expression results in a
valid type, `sfinae_true` can be instantiated and inherits from
`std::true_type`, which provides the `::value = true` static member.

Next, we will follow the instantiation chain of `hasSize` for the case in
which Container does implement a size member.

- We instantiate `hasSize<ContainerWithSize>`. This struct tries to inherit
  from the type specified by `detail::testHasSize<ContainerWithSize>(0)`. There
  are two functions of the required name where the literal `0` is implicitly
  convertible to the function argument: the one with an `int`, and the one
  with a `long` argument.
- The `int` version goes first, since it is declared first. It, in turn,
  attempts to inherit from the instantiation of `sfinae_true` with the template
  argument that is the type resulting from the expression
  `std::declval<ContainerWithSize>().size()`. Since our `ContainerWithSize`
  type does implement size, this expression does not result in an expression
  failure, but in the type that our `size()` member returns. Hence,
  `sfinae_true` can be instantiated, inherits from `std::true_type`, and
  following the inheritance chain, `hasSize<ContainerWithSize>::value = true`.

In the alternate instantiation chain where `Container` does not implement a
`size()` member, the integer test's template argument to the instantiation of
`sfinae_true` will result in an expression failure. The integer function is then
excluded from substitution, but does not result in a compilation error yet,
since the `long` version may still be viable (SFINAE). Since the `long` version
unconditionally inherits from `std::false_type`, this works, and
`hasSize<ContainerWithoutSize>::value = false`.
