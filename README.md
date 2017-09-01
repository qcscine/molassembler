# Symmetry definitions and properties library
## Overview

Contains a number of idealized symmetry definitions for the classification of
molecular geometries. Computations can either be performed on constexpr objects
at compile-time or a collected dynamic container at runtime.

## Data provided per Symmetries

- Number of binding sites
- String representation
- Angle function returning idealized angles between symmetry positions
- Idealized coordinates
- A minimal set of rotations to allow the generation of all superimposable
  sequences of binding site mapping indices
- A set of tetrahedra whose signed volumes permit the distinction of chiral
  elements within the symmetry


## Properties / Algorithms

Properties are implemented in both a dynamic and a constexpr fashion to ensure
correctness / consistency.

- Apply rotation to a sequence of indices
- Calculate a signed tetrahedron volume in 3D cartesian space
- Find periodicity of rotations
- Generate all superimposable index sequences for a particular starting index
  sequence within a specific symmetry
- Find ideal (meaning least angular and chiral distortion) mappings between
  symmetry positions when adding or removing a ligand as well as when a
  transition between symmetries of same size are sought.
