# Template Magic library
## Overview

A collection of useful classes and functions designed for runtime application.


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

- Some stubs to work with numeric data
  - Non-compensated summation, Kahan summation
  - Average, geometric mean, standard deviation
  - Min / Max

- Some SFINAE container and function traits
  - test for insert / push back / emplace / emplace back members
  - C++14 is callable implementation (in lieu of C++17's standardization)
  - get container value type
  - get function return type

- A minor boost::optional helper function for the construction of comparison
  operators of custom data types


## Integrating

All code is header-only due to the heavily templated nature and requires C++14.

Dependencies:

- boost: optional, any, test
- ConstexprMagic: BTree properties, Container traits


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
