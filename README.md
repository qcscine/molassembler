# Molassembler library
## Overview

Molassembler is a C++ library that aims to facilitate crossings between
Cartesian and graph representations of molecules. It provides the necessary
functionality to represent a molecule as a graph, modify it in graph space, and
generate coordinates from graphs. Local geometries are modelled onto every
non-terminal atom in the graph.


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
- Ranking algorithms are nearly fully IUPAC Blue Book 2013 compliant, extended
  to larger symmetries.
- Conformer generation uses Distance Geometry

  - Unassigned stereocenters are randomly assigned from relative statistical
    occurrence weights.
  - Full metrization during distance matrix generation scales with approximately
    N^3.5. Achieved via shortest paths calculation in graph using GOR1 algorithm.
  - Can optionally choose four-atom or 10% partial metrization.
  - Embedding and refinement is performed in four spatial dimensions.

## Integrating

This library requires the C++14 standard and a recent compiler (GCC >= 7).

External library dependencies:

- Boost 1.64
- Eigen: vector arithmetic
- dlib: BFGS solver
- Delib: Common chemical exchange types


External libraries included in tree and distribution:

- RingDecomposerLib[^1]: Unique Ring Family[^2] cycle detection
- Outcome (until released in boost): Improved error propagation
- nlohmann/json: JSON serialization
- breathe: sphinx extension to extract docstrings from doxygen xml


## Compilation

To build, run these commands starting at the main directory. 

```bash
$ mkdir build-release
$ cd build-release
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=<scine-installation-root> ..
$ make
```


## Tests

We recommend running the tests in a release build of the library. The tests are
split into two sets. After building the library and tests, run:

- Brief end-to-end tests: `make check`
- Library validation: `make validate`
- Both: `make test`

## Documentation

There are two types of documentation available: A full technical documentation
generated with Doxygen, and a tutorial-like introductory documentation built
with sphinx.

You can build both with:

```bash
$ cd documentation
$ doxygen # builds the technical documentation
$ make html # builds a tutorial-like documentation using sphinx
```

Partial technical documentations can also be built from within sub-libraries'
main folders.


[^1]: Flachsenberg, F.; Andresen, N.; Rarey, M. RingDecomposerLib: An
  Open-Source implementation of Unique Ring Families and Other Cycle Bases. J.
  Chem. Inf.  Model., 2017, 57 (2), pp 122–126

[^2]: Kolodzik, A.; Urbaczek, S.; Rarey, M. Unique Ring Families: A Chemically
  Meaningful Description of Molecular Ring Topologies. J. Chem. Inf. Model.,
  2012, 52 (8), pp 2013–2021
