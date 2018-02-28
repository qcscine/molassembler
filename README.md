# Cyclic Polygons library
## Overview

Implements methods for finding the internal angles of convex cyclic polygons
for a given set of edge lengths. Triangles are trivial, quadrilaterals still
fairly easy, but anything larger is a bit more challenging. 

Inside the library, you'll find two basic methods of finding the circumradius
of the pentagon. The first is a root-finding method for an equation from a
fairly recent paper, which appealed to me due to brevity and ease of
implementation compared to other formulations of the same equation in other
papers.

Svrtan, Dragutin. "On Circumradius Equations of Cyclic Polygons", 2009

The other is a basic central-angle summation approach, which I previously used
only to verify correctness of circumradii found in Svrtan's equation. While
Svrtan's equation also contains roots for edge-crossing non-convex cyclic
polygons, these are not of interest to us, so a simpler approach is appealing.
This method is also generalizeable to larger polygons, in contrast to Svrtan's 
equation, which is explicitly only valid for pentagons.

With the circumradius of the cyclic polygons, internal angles are directly
calculable.


## Important notes on reliability

For some reason that I cannot discern, using Svrtan's equation and some custom
root-finding algorithm is not 100% reliable if the set of edge lengths have some
significant deviation. A plot of this problem is contained in the files. The
central angle deviation method is reliable for any valid composition of edge
lengths and significantly cheaper as an added bonus.

A description at length of the two methods and the reliability issue is
presented in the documentation of CyclicPolygons.h.


## Incorporation into your project

The library is available in two forms:

- Standard linking: Treats cyclic polygons of size 3-5 only, implements both
  forms of circumradius finding.
- Header-only, minimal: Treats any cyclic polygon size. Implements only the
  central angle summation method.

Full library dependencies:
- boost: optional, test
- Temple: Composability shorthands and numeric data handling
- Constable: Compile-time square-root

Minimal dependencies:
- boost: for stable numerical root-finding
- TemplateMagic: Composability shorthands and numeric data handling

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
