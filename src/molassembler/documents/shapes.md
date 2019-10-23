# Spherical shapes (all vertices have same distance to origin)

# Two-dimensional
Line
Bent
EquilateralTriangle
T
Square
Pentagon
Hexagon

# Regular polyhedra
Tetrahedron
TrigonalPrism -> Uniform trigonal prism (square and equilateral triangle faces)
SquareAntiprism
Cube
Icosahedron
Cuboctahedron

# Reductions of regular polyhedra
VacantTetrahedron (previously apical trigonal pyramid) -> from Tetrahedron
EdgeCenteredTetragonalDisphenoid / EquatorialVacantTrigonalBipyramid -> from TrigonalBipyramid
TrigonalPyramid -> from TrigonalBipyramid

# Johnson shapes (already spherical)

SquarePyramid -> J1
PentagonalPyramid -> J2
TrigonalBipyramid -> J12
PentagonalBipyramid -> J13
HexagonalBipyramid -> Not a Johnson solid, but a bipyramid continuation
HeptagonalBipyramid -> Not a Johnson solid, but a bipyramid continuation

# non-spherical Johnson-derived shapes
CappedSquareAntiprism -> s-J10, C4v
BicappedSquareAntiprism -> s-J17, D4h
TrigonalDodecahedron -> s-J84, D2d
CappedTrigonalPrism -> s-J49, C2v
BicappedTrigonalPrism -> s-J50, C2v (drop, too close to square antiprism)
TricappedTrigonalPrism -> s-J51, D3h

To get to a spherical shape, minimize the sum of inverse distances between all
pairs of vertices. Just numerically minimizing that gets you a local minimum of
the Thomson potential which should still be point group symmetric.

# Unknown
CappedOctahedron -> C3v Gyroelongated triangular pyramid / "capped triangular antiprism"
