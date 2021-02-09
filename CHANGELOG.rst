Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <http://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <http://semver.org/spec/v2.0.0.html>`_.

1.1.0 [in progress]
-------------------

Added
.....

- Molecule: Bond stereopermutator addition and removal functions
- Molecule: Manual atom stereopermutator thermalization control
- Graph: Added modifying functions to the interface
- Graph algorithm: Ranking equivalent groups and ranking distinct atoms
  algorithms exploiting ranking results to classify chemical equivalence.
- Graph algorithm: Shortest path generator between vertices in the graph
- Graph algorithm: Edit distance. Algorithm to calculate minimal set of vertex
  and edge alterations to transform one graph into another. 
- Graph algorithm: Reaction edit distance. Variation of the graph edit distance
  algorithm that conserves element types. Also adds an associated function that
  plots the edits.
- Conformer deduplication: More ``Relabeler``-related functions
- Experimental SMILES emitter: ``IO/SmilesEmitter.h``
- Python bindings:

  - Added modifying functions to ``Graph``
  - Direct copying support for ``Molecule`` instead of via pickling
  - Added build-time generation of typing stubs with pybind11-stubgen

Changed
.......

- Graph and Molecule are now implementations of a common interface for graph
  information and modification. This interface might be expanded in the future.
- More colorful graphviz rendering for elements
- Conan: Better integration with community packages. No longer require full CMake
  from all dependencies but follow conan `packaging philosophy <https://github.com/conan-io/conan-center-index/blob/master/docs/faqs.md#why-are-cmake-findconfig-files-and-pkg-config-files-not-packaged>`_.
- CMake: option ``MOLASSEMBLER_PARALLELIZE`` is now ``SCINE_PARALLELIZE`` to
  follow SCINE convention.
- Interpretations without bond orders infer bond orders in binary fashion from
  Utils' BondDetector instead of from UFF parameters
- Bond stereopermutator alignment:
  ``BondStereopermutator::Alignment::BetweenEclipsedAndStaggered`` now generates
  the same amount of alignments as
  ``BondStereopermutator::Alignment::Eclipsed``, not twice as many.
- Auxiliary library ``Temple``

  - Refactor ``map`` to be able to apply it to tuples and arrays, too.
  - Clean up ``ContainerTraits.h``

- SMILES parser

  - Parsing errors as part of exception string, not written to stdout
  - Add partial handling of aromatic atom types in valence filling and bond
    stereopermutation instation. No kekulization validation yet.
  - Fix valence filling bug for atom types with multiple valid valences

- Python bindings

  - Altered name of ``ChiralStatePreservation`` enum member from ``None`` to
    ``DoNotPreserve`` (the former is a reserved keyword)
  - Better automatic type signature annotations in docstrings

Deprecated
..........

- ``Graph`` properties ``N`` and ``B`` for the number of atoms and the number of
  bonds have been deprecated in favor of new ``V`` and ``E`` properties in order
  to match complexity annotations and single-letter object properties
- ``StereopermutatorList`` method ``try_remove`` is deprecated in favor of
  ``remove``, which now behaves as ``try_remove`` would (no throwing).

Fixed
.....

- Conformer generation: Removed an incorrect check for non-terminal vertices
  without an atom stereopermutator
- Permutator propagation: 

  - Fixed missing propagation of atom stereopermutator placement and re-keying
    the atom stereopermutator map in StereopermutatorList
  - Add missing propagation of bond stereopermutator state on vertex removal

- Directed conformer generation: Fixed incorrect precondition check with
  unassigned stereopermutators



1.0.0
-----

Added
.....

- Add Conan support
- Explicit definition of which headers make up the public API in the
  tutorial-like documentation
- Molecule canonicalization: After canonicalization, isomorphism checks reduce
  to an identity comparison.
- GraphAlgorithms.h for public graph algorithms. Currently contains only a
  graph distance BFS algorithm 
- A ``doc`` target that builds the Doxygen documentation, which is now more
  extensive and contains the beginnings of a tutorial
- Many parameters of Distance Geometry can now be altered by passing a
  non-defaulted ``Configuration`` object.
- Isomer predicate and generator header ``Isomers.h``
- Higher-level editing functionality in ``Editing.h``
- More shapes up to icosahedron and cuboctahedron (size 12)
- Continuous symmetry, shape measures
- Experimental SMILES Parser
- Python bindings

  - Molecule instances integrate nicely with notebooks using ``_repr_svg_``
  - Doctested examples

Changed
.......

- The PRNG Engine is seeded directly instead of a wrapper object that helps
  with generating random numbers. The PRNG engine is part of molassembler's
  public interface instead of the sublibrary temple. The engine is constructed
  on first use.
- Molassembler's validation and analysis binaries are no longer built by
  default (see CMake options)
- Add BSD-3 license marker to all files and a checker script
- Adopt Scine code conventions regarding namespace formatting
- Enclose temple, shapes, and stereopermutation sub-libraries in molassembler
  namespace
