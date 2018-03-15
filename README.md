# Cyclic Polygons library
## Overview

Implements methods for finding the internal angles of convex cyclic polygons
for a given set of edge lengths. Triangles are trivial, quadrilaterals still
fairly easy, but anything larger is a bit more challenging. 

General cyclic polygon circumradius calculation is done via a basic
central-angle summation approach. Since it is difficult to determine whether the
circumcenter of the convex polygon will be inside or outside the polygon ahead
of time, both cases are treated separately.


## Incorporation into your project

The library is a single header file.

Header dependencies:
- boost: for stable numerical root-finding
- temple: Composability shorthands and numeric data handling

The code requires the C++14 standard to compile.


## Compiling and running the tests

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test
```

## Documentation

You can build the documentation by running `doxygen` in the main directory.
