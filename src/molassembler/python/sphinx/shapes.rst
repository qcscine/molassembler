Shapes
======

.. figure:: images/all_shapes.png

   All shapes that molassembler can classify from positions and treat
   permutationally. These are shapes from 2 vertices (line, bent) up to 12
   vertices (icosahedron, cuboctahedron).

   Shapes are classified into molecules by determining the nearest shape with
   a matching number of vertices.

.. figure:: images/molecule_shapes.png

   Octahedron, tetrahedron and equilateral triangle shapes shown in a sample
   molecule.

   When it comes to haptic ligands, the set of haptically bonded atoms' centroid
   is taken as a putative vertex position.

.. figure:: images/haptic_shape.png

   The centroid of the benzene in is shown in a ghostly black in the center of
   the cycle. The shape of the iron center is then closest to an equilateral
   triangle.

.. automodule:: molassembler.shapes
   :members:
   :undoc-members:
