# Temple: A Utilities library
## Overview

A collection of composability improvements for terser programming, patterns
encountered multiple times, and constexpr containers and algorithms for
compile-time programming.


## Containers

- A helper class for RNG
- A vector view class that allows filtering and sorting without modifying the
  underlying data
- A range-for enumeration adapter, providing the data along with an index for
  any type of container that implements iterators
- An adapter class for fetching the member value of containers containing
  classes with referential access instead of by copy using map
- A highly customizable B-Tree
- Two Cache implementations, one minimal with an underlying map of T -> U, one
  more involved with a boost::any ValueType

## Constexpr containers (with constexpr iterators)
- Array: like std::array, fixed-size 
- DynamicArray: managed array with fixed-size maximal number of elements, more
  like std::vector but without reallocation on exceeding reserved size
- Set: fixed-size minimal set, insertion or deletion changes type signature
- DynamicSet: fixed-maximum-size managed set, insertion or deletion does not
  change type signature
- DynamicMap: fixed-maximum-size associative container
- BTree: A B-tree that stores merely values, not key-value pairs
- UIntArray, DynamicUIntArray: Array-imitating containers that store short
  sequences of the digits 0-9 as a singular number for the benefit of fast 
  comparison operators
- Vector: 3D cartesian coordinate system constexpr class with basic
  computations
- UpperTriangularMatrix: stores data of an upper triangular matrix in linear
  array
- Optional: To represent operations that may frequently return no result
- Pair: Minimal replacement for std::pair, where some member functions are not 
  yet marked constexpr


## Algorithms

- A collection of pseudo-functional-style composable container manipulation
  functions
  - map (of map-reduce fame) with specializations for std::array and std::map
  - mapSequentialPairs and mapAllPairs
  - zipMap
  - reduce / accumulate
  - count
  - variadic concatenate
  - cast
  - reverse
  - condense an iterable to a delimiter-separated string of its values
  - group values of a container by mapping value or equality
  - copyIf
  - all of, any of
  - composable set union, intersection, difference
  - ...

- Some SFINAE container and function traits
  - test for insert / push back / emplace / emplace back members
  - C++14 is callable implementation (in lieu of C++17's standardization)
  - get container value type
  - get function return type

- boost::optional composition syntactic sugar

## Constexpr algorithms

- Containers: Some functional-style computations on containers, e.g.:

  - map, reduce, sum 
  - array concatenate, comparison, push and pop
  - lower bound binary search in ordered containers
  - insert into sorted container
  - in-place swap, reverse
  - next and previous permutation

- TupleType: type-level computations on types in a tuple

  - forward types to template function
  - map types to unary template function (also mapAllPairs to binary template
    function)
  - allOf
  - count types in tuple types

- Some stubs to work with numeric data
  - Uncompensated summation
  - Average, geometric mean, standard deviation
  - Min / Max


## Integrating

All code is header-only due to its heavily templated nature and requires C++14.

Dependencies: 
- boost


## Compiling and running tests

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test
```


## Documentation

You can build the documentation by running `doxygen` in the main directory.
