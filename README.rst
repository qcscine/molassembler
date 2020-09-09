====================
SCINE - Molassembler
====================

Overview
========

Molassembler is a C++ library that aims to facilitate conversions between
Cartesian and graph representations of molecules. It provides the necessary
functionality to represent a molecule as a graph, modify it in graph space, and
generate new coordinates from graphs. It can capture the absolute configuration
of inorganic molecules with multidentate and haptic ligands from Cartesian
coordinates and enumerate non-superposable stereopermutations at non-terminal
atoms and non-isotropic bonds at arbitrary local shapes ranging up to the
icosahedron and cuboctahedron.

Features
========

- Molecules can be constructed from many types of information.
- Stereocenters are treated in shapes ranging from monovacant tetrahedron all
  the way up to the icosahedron and cuboctahedron.
- A high-temperature approximation is invoked by default to avoid considering
  inverting nitrogen centers as stereocenters, but this is optional. Even in
  the high-temperature approximation, nitrogen centers whose substituents
  form a strained cycle and hence do not invert rapidly are considered a
  stereocenter.
- All stereocenter permutations are generated with relative statistical
  occurrence weights. Linking of ligands (denticity) is properly considered.
  Several classes of haptic ligands are supported.
- Editing of molecules preserves chiral information by default, and is highly
  configurable.
- Molecules can be canonicalized for fast isomorphism tests. Canonicalization
  can be customized to use subsets of the available information for vertex
  coloring if desired.
- Ranking algorithms are nearly fully IUPAC Blue Book 2013 compliant, extended
  to larger symmetries.
- Stochastic conformer generation with Distance Geometry
- Directed conformer generation through enumeration of rotamers

License
=======

Molassembler is licensed under the BSD 3-clause "New" or "Revised" license. See
also the ``LICENSE.txt`` file.

Integrating
===========

This library requires the C++14 standard.

Dependencies:

- SCINE Utils (BSD-3 license) >= 3.0.0
- Boost (Boost license) >= 1.65 (lowest tested, prefer newest)
- Eigen (MPL 2.0 license) >= 3.3.2
- (BLAS library, added if detected during compilation)

Can currently be compiled with:

- [x] GCC >= 7
- [x] Clang >= 4
- [x] MinGW-w64 (latest)
- [ ] MSVC (compiler compliance issues with ``constexpr``)

Unowned libraries included in distribution (see ``src/extern``):

- RingDecomposerLib [1]_ (BSD-3 license): Unique Ring Family [2]_ cycle detection
- Outcome (until released in boost): Improved error propagation
- nlohmann/json (MIT license): JSON serialization
- nauty [3]_ (Apache 2.0 license): Graph automorphism determination and canonical labeling

This library uses CMake to model dependencies and make builds
platform-independent.

How to Cite
===========

When publishing results obtained with Molassembler, please cite the
corresponding release as archived on Zenodo **TODO**.

In addition, we kindly request you cite the following article when using
Molassembler:

J.-G. Sobez, M. Reiher, "Molassembler: Molecular Graph Construction,
Modification, and Conformer Generation for Inorganic and Organic
Molecules", *J. Chem. Inf. Model*, **2020**, *60*, 3884.

Compilation
===========

To build with CMake only, run these commands starting at the main directory::

    mkdir build-release
    cd build-release
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

To build with Conan::

    conan create --build=missing .

Tests
=====

We recommend running the tests in a release build of the library. The debug
builds can run for a good 10 minutes. After building the library and tests,
run ``make test``. The Python bindings are tested with ``pytest`` and ``doctest``,
if available.


Documentation
=============

Built documentation for releases is available for the `C++ library`_ and `Python bindings`_.

If ``doxygen`` is found, the C++ library documentation is built. If the Python
bindings are built and the ``sphinx`` Python module is available, the Python
binding documentation is generated too.

.. _C++ library: https://scine.ethz.ch/static/download/documentation/molassembler/v1.0.0/cpp/index.html

.. _Python bindings: https://scine.ethz.ch/static/download/documentation/molassembler/v1.0.0/py/index.html

.. [1] Flachsenberg, F.; Andresen, N.; Rarey, M. RingDecomposerLib: An
       Open-Source implementation of Unique Ring Families and Other Cycle Bases. J.
       Chem. Inf. Model., 2017, 57 (2), pp 122–126

.. [2] Kolodzik, A.; Urbaczek, S.; Rarey, M. Unique Ring Families: A Chemically
       Meaningful Description of Molecular Ring Topologies. J. Chem. Inf. Model.,
       2012, 52 (8), pp 2013–2021

.. [3] McKay, Brendan D.; Adolfo Piperno. Practical graph isomorphism, II.
       J. Symb. Comput., 2014, 60, pp 94-112.
