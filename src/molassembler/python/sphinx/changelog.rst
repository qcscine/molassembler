Changelog
=========

This project adheres to semantic versioning.

1.0.0 - 2019-12-01
------------------
Added
~~~~~
- Explicit definition of which headers make up the public API in the
  tutorial-like documentation
- Molecule canonicalization: After canonicalization, isomorphism checks reduce
  to an identity comparison.
- GraphAlgorithms.h for public graph algorithms. Currently contains only a
  graph distance BFS algorithm 
- A `doc` target that builds the Doxygen documentation, which is now more
  extensive and contains the beginnings of a tutorial
- Many parameters of Distance Geometry can now be altered by passing a
  non-defaulted `Configuration` object.
- Isomer predicate and generator header `Isomers.h`
- Higher-level editing functionality in `Editing.h`
- Python bindings
  - Molecule instances integrate nicely with notebooks using `_repr_svg_`

Changed
~~~~~~~
- temple::Bitmask no longer forms part of the public interface
- Rebased the library on Scine's UtilsOS instead of on the Delib (which is now
  deprecated)
- Renamed AtomStereocenter and BondStereocenter to AtomStereopermutator and
  BondStereopermutator respectively. StereocenterList is accordingly renamed to
  StereopermutatorList. The classes in question often handle cases where there
  is merely a single stereopermutation and hence **are not** a stereocenter. It
  seemed misleading to keep the name. 
- Removed setModelInformation and chiralityConstraints from Stereopermutator
  interfaces. This is now in sole custody of SpatialModel, and should not affect
  any consumers.
- The PRNG Engine is seeded directly instead of a wrapper object that helps
  with generating random numbers. The PRNG engine is part of molassembler's
  public interface instead of the sublibrary temple. The engine is constructed
  on first use.
- Improved CMake code to better match modern use (remove `GLOB`, do not pollute
  `CMAKE_CXX_FLAGS`, etc.)
- Molassembler's validation and analysis are no longer built by default (see
  CMake options)
- Library is no longer based on Delib exchange formats, but rather on scine
  utils open source library (which are virtually identical).
