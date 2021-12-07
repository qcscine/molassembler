Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <http://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <http://semver.org/spec/v2.0.0.html>`_.

1.1.0 [in progress]
-------------------

Added
.....

- Graph algorithm: Ranking equivalent groups and ranking distinct atoms
  algorithms exploiting ranking results to classify chemical equivalence.
- Graph algorithm: Shortest path generator between vertices in the graph
- Conformer deduplication: More ``Relabeler``-related functions
- Python bindings:

  - Added modifying functions to ``Graph``
  - Direct copying support for ``Molecule`` (``__copy__`` and ``__deepcopy__``)
    instead of indirect via pickling

Changed
.......

- Conan: Better integration with community packages. No longer require full CMake
  from all dependencies but follow conan `packaging philosophy <https://github.com/conan-io/conan-center-index/blob/master/docs/faqs.md#why-are-cmake-findconfig-files-and-pkg-config-files-not-packaged>`_.
- CMake: option ``MOLASSEMBLER_PARALLELIZE`` is now ``SCINE_PARALLELIZE`` to
  follow SCINE convention.
- Bond stereopermutator alignment:
  ``BondStereopermutator::Alignment::BetweenEclipsedAndStaggered`` now generates
  the same amount of alignments as
  ``BondStereopermutator::Alignment::Eclipsed``, not twice as many.
- Auxiliary library ``Temple``

  - Refactor ``map`` to be able to apply it to tuples and arrays, too.
  - Clean up ``ContainerTraits.h``

Deprecated
..........

Removed
.......

Fixed
.....

- Conformer generation: Removed an incorrect check for non-terminal vertices
  without an atom stereopermutator
- Permutator propagation: 

  - Fixed missing propagation of atom stereopermutator placement and re-keying
    the atom stereopermutator map in StereopermutatorList
  - Add missing propagation of bond stereopermutator state on vertex removal


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
