# Constexpr Utilities library
## Overview

Contains a number of useful constexpr containers and algorithms for compile-time
programming.


## Containers

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
