Stereopermutators
=================

Currently, this library models two types of stereocenters: The first models the
atom-centric relative placement of substituents across geometries ranging from
linear to square-antiprismatic. The second type models configurational
permutations owing to bond-centric rotational barriers.

In order to determine whether a particular atom in a molecule is an atom-centric
stereocenter, its substituents are ranked according to a nearly-compliant
implementation of the IUPAC sequence rules as laid out in the 2013 Blue Book.
The rankings are transformed into an abstract substituent case (e.g. octahedral
(A-A)BBCD) and a symbolic computation is carried out to determine the number of
permutations that are not superimposable via spatial rotations within the
idealized shape. The set of resulting permutations of the substituent symbols
is called the set of stereopermutations. If this set contains more than one
stereopermutation, then the atom is an atom-centric stereocenter under that
shape.

If substituents are haptic or multidentate, an additional algorithm removes
stereopermutations it deems clearly impossible. All bridge lengths between
pairs of chelating atoms of a multidentate ligand are checked against the atom
pair's bite angle within the idealized shape. Additionally, haptic ligands'
cones are checked to ensure they do not overlap. Indices within the set of not
clearly impossible stereopermutations are called assignments.

Bond-centric stereocenters model rotational configurations of arbitrary
combinations of two shapes and their respective fused shape positions.
The fused shape positions of each side affect the overall permutations if the
shape has multiple position groups. For instance, this is the case in square
pyramidal shapes, where there are axial and equatorial shape positions.

.. autoclass:: scine_molassembler.StereopermutatorList
   :members:
   :undoc-members:

Atom stereopermutators
----------------------

.. autoclass:: scine_molassembler.AtomStereopermutator
   :members:
   :undoc-members:

Bond stereopermutators
----------------------

.. autoclass:: scine_molassembler.BondStereopermutator
   :members:
   :undoc-members:
