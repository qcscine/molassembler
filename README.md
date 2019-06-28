# Molassembler library
## Overview

Molassembler is a C++ library that aims to facilitate crossings between
Cartesian and graph representations of molecules. It provides the necessary
functionality to represent a molecule as a graph, modify it in graph space, and
generate coordinates from graphs. It can capture the absolute configuration
of multidentate and haptic inorganic molecules from positional data and
generate non-superposable stereopermutations as output.


## Features

- Molecule construction from many types of information. 
- Stereocenters are treated from trigonal pyramidal all the way up to square
  antiprismatic local geometries.
- High-temperature approximation is invoked by default to avoid considering
  inverting nitrogen centers as stereocenters, but this is optional. Even in
  the high-temperature approximation, nitrogen centers whose substituents
  form a strained cycle and hence do not invert rapidly are considered a
  stereocenter.
- All stereocenter permutations are generated with relative statistical
  occurrence weights. Linking of ligands (denticity) is properly considered.
  Several classes of haptic ligands are supported.
- Editing of molecules preservers chiral information by default, and is highly
  configurable.
- Molecules can be canonicalized for fast isomorphism tests. Canonicalization
  can be customized to use subsets of the available information for vertex
  coloring if desired.
- Ranking algorithms are nearly fully IUPAC Blue Book 2013 compliant, extended
  to larger symmetries.
- Stochastic conformer generation with Distance Geometry
  - Unassigned stereocenters are randomly assigned from relative statistical
    occurrence weights.
  - Full metrization during distance matrix generation scales with approximately
    N^3.5. Achieved via shortest paths calculation in graph using GOR1 algorithm.
  - Can optionally choose four-atom or 10% partial metrization.
  - Embedding and refinement is performed in four spatial dimensions.


## Integrating

This library requires the C++14 standard.

Dependencies:

- Boost >= 1.64
- Eigen >= 3
- SCINE Utils > 0.1.0
- (BLAS library, added if detected during compilation)


Can currently be compiled with:

- [x] GCC >= 7
- [x] Clang >= 4
- [ ] MSVC (compiler compliance issues with C++14 `constexpr` code)

Windows compatibility is in progress. In the meantime, consider options like
MinGW (compiles with GCC) or Visual Studio Codegen with Clang to create Windows
libraries.

Unowned libraries included this distribution:

- RingDecomposerLib[^1]: Unique Ring Family[^2] cycle detection
- Outcome (until released in boost): Improved error propagation
- nlohmann/json: JSON serialization
- spectra: ARPACK-like interface built on top of Eigen
- nauty[^3]: Graph automorphism determination and canonical labeling

This library uses CMake to model dependencies and make builds
platform-independent.


## Compilation

To build, run these commands starting at the main directory. 

```bash
$ mkdir build-release
$ cd build-release
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```

## Tests

We recommend running the tests in a release build of the library. The tests are
split into two sets. After building the library and tests, run:

- Brief end-to-end tests: `make check`
- Library validation: `make validate`
- Both: `make test`


## Documentation

The documentation can be built from the root directory with:

```bash
$ cd build
$ cmake ..
$ make doc
```

This generates the Doxygen documentation in a `doc` subfolder.

[^1]: Flachsenberg, F.; Andresen, N.; Rarey, M. RingDecomposerLib: An
  Open-Source implementation of Unique Ring Families and Other Cycle Bases. J.
  Chem. Inf.  Model., 2017, 57 (2), pp 122–126

[^2]: Kolodzik, A.; Urbaczek, S.; Rarey, M. Unique Ring Families: A Chemically
  Meaningful Description of Molecular Ring Topologies. J. Chem. Inf. Model.,
  2012, 52 (8), pp 2013–2021

[^3]: McKay, Brendan D.; Adolfo Piperno. Practical graph isomorphism, II.
  J. Symb. Comput., 2014, 60, pp 94-112.
