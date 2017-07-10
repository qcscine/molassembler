# Symmetry information library
## Overview

Provides structured information about a set of symmetries. For every symmetry, 
the following information is exposed:

- A string name
- The number of substituents / coordination number
- A minimal set of rotations that allows the construction of any rotational
  equivalent
- A function that provides the geometric angle between any two unequal
  substituents
- A list of tetrahedron indices that can encapsulate all steric information 
  contained in the symmetry
- A set of sample coordinates that fulfill the set geometric criteria

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

## Integrating

Is minimally dependent on boost and Eigen. If you choose to use \c constexpr
precalculated angles for the square antiprismatic symmetry for better angle
function correctness, the library is also dependent on ConstexprMagic.
