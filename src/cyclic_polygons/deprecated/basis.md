
Inside the library, you'll find two basic methods of finding the circumradius
of the pentagon. The first is a root-finding method for an equation from a
fairly recent paper, which appealed to me due to brevity and ease of
implementation compared to other formulations of the same equation in other
papers.

Svrtan, Dragutin. "On Circumradius Equations of Cyclic Polygons", 2009

## Important notes on reliability

For some reason that I cannot discern, using Svrtan's equation and some custom
root-finding algorithm is not 100% reliable if the set of edge lengths have some
significant deviation. A plot of this problem is contained in the files. The
central angle deviation method is reliable for any valid composition of edge
lengths and significantly cheaper as an added bonus.

A description at length of the two methods and the reliability issue is
presented in the documentation of CyclicPolygons.h.


