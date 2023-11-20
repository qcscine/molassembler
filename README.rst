SCINE - Molassembler
====================

Overview
--------

Molassembler is a C++ library that aims to facilitate conversions between
Cartesian and graph representations of molecules. It provides the necessary
functionality to represent a molecule as a graph, modify it in graph space, and
generate new coordinates from graphs. It can capture the absolute configuration
of inorganic molecules with multidentate and haptic ligands from Cartesian
coordinates and enumerate non-superposable stereopermutations at non-terminal
atoms and non-isotropic bonds at arbitrary local shapes ranging up to the
icosahedron and cuboctahedron.

Features
--------

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
- Ranking algorithms are nearly fully IUPAC Blue Book 2013 compliant,
  generalized to larger coordination polyhedra.
- Stochastic conformer generation with Distance Geometry
- Directed conformer generation through enumeration of rotamers

License
-------

Molassembler is licensed under the BSD 3-clause "New" or "Revised" license. See
also the ``LICENSE.txt`` file.

Integrating
-----------

This library requires the C++17 standard.

Dependencies:

- `SCINE Utilities <https://github.com/qcscine/utilities>`_ (BSD-3 license) >= 3.0.0
- `Boost <https://www.boost.org/>`_ (Boost license) >= 1.65 (lowest tested, prefer newest)
- `Eigen <http://eigen.tuxfamily.org>`_ (MPL 2.0 license) >= 3.3.2
- `RingDecomposerLib <https://github.com/rareylab/RingDecomposerLib>`_ [1]_ (BSD-3 license): Unique Ring Family [2]_ cycle detection
- `Outcome <https://github.com/ned14/outcome>`_ single-header (Boost license): Enforce error handling requirement in type system
- `JSON For Modern C++ <https://github.com/nlohmann/json>`_ (MIT license): JSON serialization
- `nauty <http://pallini.di.uniroma1.it>`_ [3]_ (Apache 2.0 license): Graph automorphism determination and canonical labeling
- (MKL/LAPACK/BLAS libraries, added if detected during compilation)

Can currently be compiled with:

- GCC >= 7
- Clang >= 4
- MinGW-w64 (latest)

MSVC is currently untested. Last attempts failed because of compiler standard
compliance issues with ``constexpr`` code.

How to Cite
-----------

When publishing results obtained with Molassembler, please cite the
corresponding release as archived on `Zenodo <https://doi.org/10.5281/zenodo.4293554>`_.

In addition, we kindly request you cite the following article when using
Molassembler:

J.-G. Sobez, M. Reiher, "Molassembler: Molecular Graph Construction,
Modification, and Conformer Generation for Inorganic and Organic
Molecules", *J. Chem. Inf. Model*, **2020**, *60*, 3884.

Installation
------------

CMake
.....

When building with CMake, Boost and Eigen must be installed and available via
CMake's ``find_package`` (e.g. via ``CMAKE_PREFIX_PATH``). Any of the other
libraries can be available, but are downloaded dynamically if missing. 

Clone the repository, then enter the following commands::

    git submodule update --init
    mkdir build-release
    cd build-release
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

You may want to peruse the CMake options to disable building the tests or
activating the Python binding builds. Run ``cmake -L ..`` to list options
affecting the build. Look for options with the ``SCINE_`` prefix.

Conan
.....

No dependencies must be preinstalled, and you do not need to download the
sources. To install/build with Conan::

    conan remote add scine https://scine-artifactory.ethz.ch/artifactory/api/conan/public
    conan install -r scine --build=missing scine_molassembler/2.0.1@

Should you want Python bindings, add ``-o scine_molassembler:python=True`` before
the last argument.

PyPI
....

``manylinux`` packages of thie Python bindings are available from PyPI and can be
installed with::

    python3 -m pip install scine_molassembler

Documentation
-------------

Built documentation for releases is available for the `C++ library`_ and `Python bindings`_.

If ``doxygen`` is found, the C++ library documentation is built. If the Python
bindings are built and the ``sphinx`` Python module is available, the Python
binding documentation is generated too.

.. _C++ library: https://scine.ethz.ch/static/download/documentation/molassembler/v2.0.1/cpp/index.html

.. _Python bindings: https://scine.ethz.ch/static/download/documentation/molassembler/v2.0.1/py/index.html

References
----------

.. [1] Flachsenberg, F.; Andresen, N.; Rarey, M. RingDecomposerLib: An
       Open-Source implementation of Unique Ring Families and Other Cycle Bases. J.
       Chem. Inf. Model., 2017, 57 (2), pp 122–126

.. [2] Kolodzik, A.; Urbaczek, S.; Rarey, M. Unique Ring Families: A Chemically
       Meaningful Description of Molecular Ring Topologies. J. Chem. Inf. Model.,
       2012, 52 (8), pp 2013–2021

.. [3] McKay, Brendan D.; Adolfo Piperno. Practical graph isomorphism, II.
       J. Symb. Comput., 2014, 60, pp 94-112.
