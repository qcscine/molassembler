# Shape definitions and properties library
## Overview

Contains a number of idealized shape definitions for the classification of
local molecular geometries. Computations can either be performed on constexpr
objects at compile-time or their dynamic counterparts at runtime. Due to the
constexpr nature of some algorithms, C++14 is required.


## Data provided per shape

- Number of binding sites
- String representation
- Angle function returning idealized angles between shape vertex positions
- Idealized coordinates
- A minimal set of rotations to allow the generation of all superimposable
  sequences of binding site mapping indices
- A set of tetrahedra whose signed volumes permit the distinction of chiral
  elements within the shape
- A mirror permutation


## Properties / Algorithms

Properties are implemented in both a dynamic and a constexpr fashion to ensure
correctness / consistency.

- Apply rotation to a sequence of indices
- Calculate a signed tetrahedron volume in 3D cartesian space
- Find periodicity of rotations
- Generate all superimposable index sequences for a particular starting index
  sequence within a specific shape
- Find ideal (meaning least angular and chiral distortion) mappings between
  shape positions when adding or removing a ligand as well as when a
  transition between shapes of same size are sought.
- ...


## Compilation

Since `constexpr` is a relatively new C++ feature, early compiler
implementations have been spotty. It is known that gcc 6.3.3 and 7.3.0 are
unable to handle the precomputation of several constexpr algorithms, whose
results would otherwise be available at runtime. Newer versions of gcc can, but
require a significant amount of time and memory.

Clang >= 4.0.0 is known to compile these algorithms just fine in 3-5 minutes.

The property precomputation is currently only enabled for the clang family of
compilers. If you nevertheless want to try the precomputation with your
particular compiler, specify `cmake -DSHAPES_TRY_CONSTEXPR ..` in your build
directory.
