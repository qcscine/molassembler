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


## Container adaptor argument conditional ownership

If you want to create a new range from some container, the handling of whether
the container being modified is owned by the new range is very important to
avoid both extraneous copies and dangling references.

We have to differentiate between rvalue and lvalue reference arguments to the
adaptor constructor. Any rvalues have to be owned by the range being
constructed, while any lvalues should be bound to a (const) reference. All
ranges are template classes, and we shall call one such template argument T. We
can compute the necessary type by applying the following transformation:

```c++ 
using BoundT = std::conditional_t<
  std::is_rvalue_reference<T&&>::value,
  std::decay_t<T>,
  const T&
>;
```

Since T is used in a type-deducing context, `T&&` applies a universal reference,
applying reference collapsing to any reference members of T:

```c++
T = std::vector<unsigned>&
T&& == std::vector<unsigned>&

T = std::vector<unsigned>&&
T&& == std::vector<unsigned>&&
```

If the result of `T&&` is an rvalue, we have to **own** the T, so we apply
`std::decay_t<T>`, which removes any const, volatile or reference type
qualifiers. If `T&&` is not an rvalue, we can bind it to a const reference.
Remember here that `T& + & == T&`.

When defining the constructor of the adaptor, we can apply the same logic
unconditionally of whether we have to own the T or not:

```c++
explicit Adaptor(Container&& passContainer)
  : boundContainer(std::forward<Container>(passContainer)) {}
```
