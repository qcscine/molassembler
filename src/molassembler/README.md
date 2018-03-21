# Molassembler library
## Overview

Molassembler is a C++ library that aims to facilitate crossings between
Cartesian and graph representations of chemical molecules. It provides the
necessary functionality to represent a chemical molecule as a graph, modify it
in graph space, and generate coordinates from graphs. Local geometries are
modelled onto every non-terminal atom in the graph.


## Features

- MOLFile input and output, stereocenter assignment detection from 3D
  coordinates if present
- Stereocenters are considered from trigonal pyramidal all the way up to square
  antiprismatic local geometries.

  - High-temperature approximation is invoked by default to avoid considering
    inverting nitrogen centers as stereocenters, but can be turned off. Even in
    the high-temperature approximation, nitrogen centers whose substituents
    form a strained cycle and hence do not invert quickly are considered a
    stereocenter.
  - All non-superimposable configurations of a stereocenter are generated with
    relative statistical occurrence weights, internally called "assignments".
    Linking of ligands (denticity) is properly considered.
  - Ligand additions, ligand removals and symmetry alterations on stereocenters
    can preserve steric information in several variants:
    - Do not preserve steric information
    - Effortless and unique mappings (default)
    - Unique mappings
    - Random from multiple best

- Ranking algorithms are nearly fully IUPAC compliant.
- Cycle detection using Unique Ring Families
- Distance Geometry algorithm to generate Cartesian coordinates

  - Unassigned stereocenters are randomly assigned from relative statistical
    occurrence weights.
  - Full metrization during distance matrix generation scales with approximately
    N^3.5. Achieved via shortest paths calculation in graph using GOR algorithm.
  - Can optionally choose four-atom or 10% partial metrization.
  - Embedding and refinement is performed in four spatial dimensions.

- Geometry determination algorithms

  - Currently only VSEPR, fallback is random symmetry of correct size

- Extensive tests

  - The ranking algorithm is tested against examples from the IUPAC
    Blue Book, the corresponding MOLFiles can be found in
    tests/mol_files/ranking_tree_molecules/.
  - The conformational model is tested against a battery of strained molecules
    which test the adaptability of geometrical idealization strictness. Those
    files can be found in tests/mol_files/strained_organic_molecules/.
  - Sanity tests for every stereocenter symmetry and assignment to ensure
    stereocenters are generated as desired
  - Various other tests.

- Analysis binaries

  - RaytraceRefinement creates POV-Ray files and .csv files with which the
    refinement stage of DistanceGeometry can be examined in detail.
    Complementary files for graphing the output and generating the ray-traced 
    output can be found in analysis/raytrace/.
  - RankingExplainer generates Graphviz files that shows individual ranking
    algorithm steps (NOTE: works only with debug build).
  - Gor1Explainer shows the progression of the simplified GOR1 single source
    shortest paths algorithm through the specially generated graph for
    DistanceGeometry triangle inequality bounds determination.
  - BenchmarkGraphAlgorithms generates benchmark data for various combinations
    of individual distances matrix generation algorithms with data types and
    metrization options. Complementary files for graph output are found in
    analysis/graph_benchmarks/


## Integrating

This library requires the C++14 standard.

External library dependencies:

- boost: graph, functional, optional, test
- Eigen: vector arithmetic
- dlib: BFGS solver
- RingDecomposerLib[^1]: Unique Ring Family[^2] cycle detection
- Delib: Common chemical exchange types

Sub-libraries:

- constable: constexpr algorithms and data structures
- temple: randomness, cache, composability improvement shorthands
- Symmetry information: constexpr properties of symmetries
- Unique Assignment: Stereopermutation of ligands in arbitrary symmetries
- Cyclic polygons: Determination of planar cyclic polygon internal angles for
  any edge lengths


## Compilation

To build, run these commands starting at the main directory. 

```bash
$ mkdir build-release
$ cd build-release
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
```


## Tests

We recommend running the tests in a debug build of the library. Since the tests
are fairly extensive, they may require considerable time to complete.

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test
```

## Documentation

You can build the documentation by running `doxygen` in the main directory.

## References

[^1]: Flachsenberg, F.; Andresen, N.; Rarey, M. RingDecomposerLib: An
  Open-Source implementation of Unique Ring Families and Other Cycle Bases. J.
  Chem. Inf.  Model., 2017, 57 (2), pp 122–126

[^2]: Kolodzik, A.; Urbaczek, S.; Rarey, M. Unique Ring Families: A Chemically
  Meaningful Description of Molecular Ring Topologies. J. Chem. Inf. Model.,
  2012, 52 (8), pp 2013–2021
