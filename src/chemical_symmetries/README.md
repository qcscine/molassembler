# Symmetry definitions and properties library
## Overview

Contains a number of idealized symmetry definitions for the classification of
molecular geometries. Computations can either be performed on constexpr objects
at compile-time or a collected dynamic container at runtime. Due to the
constexpr nature of the algorithms, C++14 is required.


## Data provided per symmetry

- Number of binding sites
- String representation
- Angle function returning idealized angles between symmetry positions
- Idealized coordinates
- A minimal set of rotations to allow the generation of all superimposable
  sequences of binding site mapping indices
- A set of tetrahedra whose signed volumes permit the distinction of chiral
  elements within the symmetry


## Contained symmetries

- Linear (2)
- Bent (2)
- Trigonal planar (3)
- Trigonal pyramidal (3)
- T-shaped (3)
- Tetrahedral (4)
- Square planar (4)
- Seesaw (4)
- Square pyramidal (5)
- Trigonal bipyramidal (5)
- Pentagonal planar (5)
- Octahedral (6)
- Trigonal prismatic (6)
- Pentagonal pyramidal (6)
- Pentagonal bipyramidal (7)
- Square antiprismatic (8)


## Properties / Algorithms

Properties are implemented in both a dynamic and a constexpr fashion to ensure
correctness / consistency.

- Apply rotation to a sequence of indices
- Calculate a signed tetrahedron volume in 3D cartesian space
- Find periodicity of rotations
- Generate all superimposable index sequences for a particular starting index
  sequence within a specific symmetry
- Find ideal (meaning least angular and chiral distortion) mappings between
  symmetry positions when adding or removing a ligand as well as when a
  transition between symmetries of same size are sought.


## Compile-time options

The behavior of the library can be changed in slight ways using compile-time
options defined in CompileTimeOptions.h. You can toggle:

- The angle function of the square antiprismatic symmetry can be pre-calculated
  as a look-up table from reference coordinates instead of from an idealized 
  representation (Default on)
- The pre-calculation of least-distortion mappings between symmetries when 
  ligands are added or removed. This can cost several (5-10) minutes at
  compile-time, but then the results are available at O(1) at run-time.
  (Default on)
- An experimental reduced set of tetrahedra for each symmetry (Default off)
- The maximum number of rotationally non-superimposable assignments for any
  symmetry as a function of how many ligands are identical, where all ligands
  are unlinked


## Integrating

This library requires the C++14 standard.

Library dependencies:

- STL
- boost: optional type, testing
- Eigen: vector arithmetic
- temple:
  - Cache
  - Composability improvement shorthands
  - constexpr algorithms and data structures


## Compilation

Since constexpr is a relatively new and complex C++ feature, compiler
implementations have been spotty. It is known that gcc 6.3.3 and 7.3.0 are
unable to handle the precomputation of several constexpr algorithms, whose
results would otherwise be available at runtime. Newer versions of gcc have not
been tested.

Clang >= 4.0.0 is known to compile these algorithms just fine in 3-5 minutes.

Since no other compilers have been tested so far, the symmetry property
precomputation is enabled only on clang builds.

To build, run these commands starting at the main directory.

```bash
$ mkdir build
$ cd build
$ export CC = path/to/clang
$ export CXX = path/to/clang++
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
$ make test
```

## Documentation

You can build the documentation by running `doxygen` in the main directory.
