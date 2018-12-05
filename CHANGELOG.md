# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

TODO
- date
- Bug fix regarding chiral constraints

## [1.0.0] - 2018-12-DD
### Added
- Explicit definition of which headers make up the public API in the
  tutorial-like documentation
- GraphAlgorithms.h for public graph algorithms. Currently contains only a
  graph distance BFS algorithm 
- A `doc` target that builds the Doxygen documentation, which is now more
  extensive and contains the beginnings of a tutorial

### Changed
- Symmetries are no longer excluded on principle, but using tau criteria (see
  AtomStereopermutator fitting)
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
- Add ETH license to all files and a checker script

### Fixed
- Several ranking bugs
- `molassembler` now correctly links against the shared `Delib` library
- `molassemblerStatic` no longer has interface dependencies on header-only
  libraries used in its implementation only
- Only permit downloading of `Delib` from GitLab if explicitly enabled for CI
- Application of IUPAC Sequence rule 5 was improved

### Removed
- The sphinx parallel documentation, including the distributed sphinx breathe
  extension to get C++ docstrings from doxygen

## [Unreleased]
### Added
### Changed
### Deprecated
### Removed
### Fixed

## [X.Y.Z] - 20YY-MM-DD
