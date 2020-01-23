# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2020-01-23
### Added
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
- More shapes up to icosahedron and cuboctahedron (size 12)
- Continuous symmetry, shape measures
- Experimental SMILES Parser
- Python bindings
  - Molecule instances integrate nicely with notebooks using `_repr_svg_`
  - Doctested examples

### Changed
- `temple` no longer forms part of the public interface
- Renamed `AtomStereocenter` and `BondStereocenter` to `AtomStereopermutator`
  and `BondStereopermutator` respectively. `StereocenterList` is accordingly
  renamed to `StereopermutatorList`. The classes in question often handle cases
  where there is merely a single stereopermutation and hence **are not** a
  stereocenter. It seemed misleading to keep the name. 
- The PRNG Engine is seeded directly instead of a wrapper object that helps
  with generating random numbers. The PRNG engine is part of molassembler's
  public interface instead of the sublibrary temple. The engine is constructed
  on first use.
- Molassembler's validation and analysis binaries are no longer built by
  default (see CMake options)
- Add BSD-3 license marker to all files and a checker script

## [Unreleased]
### Added
### Changed
### Deprecated
### Removed
### Fixed

## [X.Y.Z] - 20YY-MM-DD
