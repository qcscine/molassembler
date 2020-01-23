molassembler
============

Molassembler is a C++ library that aims to facilitate conversions between
Cartesian and graph representations of molecules. It provides the necessary
functionality to represent a molecule as a graph, modify it in graph space, and
generate new coordinates from graphs. It can capture the absolute configuration
of inorganic molecules with multidentate and haptic ligands from Cartesian
coordinates and enumerate non-superposable stereopermutations at non-terminal
atoms and non-isotropic bonds at arbitrary local shapes ranging up to the
icosahedron and cuboctahedron.

This is the documentation for molassembler's Python bindings which bind the core
functionality of the library for quick prototyping.

Core features
-------------

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
- Editing of molecules preservers chiral information by default, and is highly
  configurable.
- Molecules can be canonicalized for fast isomorphism tests. Canonicalization
  can be customized to use subsets of the available information for vertex
  coloring if desired.
- Ranking algorithms are nearly fully IUPAC Blue Book 2013 compliant, extended
  to larger symmetries.
- Stochastic conformer generation with Distance Geometry
- Directed conformer generation through enumeration of rotamers

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Introduction

   model

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: API Documentation

   molassembler
   molecule
   graph
   cycles
   stereopermutators

   conformers
   directed_conformer_generator
   editing
   interpret
   io
   options
   ranking_information
   serialization
   shapes

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Other

   issues
   changelog
   C++ API <https://scine.ethz.ch/download/molassembler/doc-cpp/index.html>
   genindex
