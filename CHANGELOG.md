# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

TODO
- date
- Bug fix regarding chiral constraints

## [1.0.0] - 2018-09-DD
### Added
- Explicit definition of which headers make up the public API in the
  tutorial-like documentation
- GraphAlgorithms.h for public graph algorithms. Currently contains only a
  graph distance BFS algorithm 

### Changed
- PRNG seeding has been changed. The PRNG Engine is seeded directly instead of
  a wrapper object that helps with generating random numbers. The PRNG engine
  is part of molassembler's public interface instead of the sublibrary temple.
  The engine is constructed on first use.
- Improved CMake code to better match modern use (remove `GLOB`, do not pollute
  `CMAKE_CXX_FLAGS`, etc.)
- Only the public API headers for `molassembler` are now installed
- Molassembler's validation and analysis are no longer built by default (see
  CMake options)
- Add ETH license to all files and a checker script

### Fixed
- `molassembler` now correctly links against the shared `Delib` library
- `molassemblerStatic` no longer has interface dependencies on header-only
  libraries used in its implementation only
- Only permit downloading of `Delib` from GitLab if explicitly enabled for CI

## [Unreleased]
### Added
### Changed
### Deprecated
### Removed
### Fixed

## [X.Y.Z] - 20YY-MM-DD
