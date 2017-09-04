# Symmetry definitions and properties library
## Overview

Contains a number of idealized symmetry definitions for the classification of
molecular geometries. Computations can either be performed on constexpr objects
at compile-time or a collected dynamic container at runtime.

## Data provided per symmetry

- Number of binding sites
- String representation
- Angle function returning idealized angles between symmetry positions
- Idealized coordinates
- A minimal set of rotations to allow the generation of all superimposable
  sequences of binding site mapping indices
- A set of tetrahedra whose signed volumes permit the distinction of chiral
  elements within the symmetry

## Contained symmetries

- Linear (2)
- Bent (2)
- Trigonal planar (3)
- Trigonal pyramidal (3)
- T-shaped (3)
- Tetrahedral (4)
- Square planar (4)
- Seesaw (4)
- Square pyramidal (5)
- Trigonal bipyramidal (5)
- Pentagonal planar (5)
- Octahedral (6)
- Trigonal prismatic (6)
- Pentagonal pyramidal (6)
- Pentagonal bipyramidal (7)
- Square antiprismatic (8)


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

## Integrating

Is minimally dependent on boost and Eigen. If you choose to use \c constexpr
precalculated angles for the square antiprismatic symmetry for better angle
function correctness, the library is also dependent on ConstexprMagic.
