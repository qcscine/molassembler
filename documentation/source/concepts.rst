========
Concepts
========

For productive use of this library, the model of various chemical concepts must
be laid plain.


Symmetries
==========

Information on the various local symmetries that are modeled are centralized
in a separate library. We have chosen to consider the following symmetries:

==== ====
Size Name
==== ====
2    Linear, bent
3    Trigonal planar, trigonal pyramidal, T-shaped
4    Tetrahedral, square planar, seesaw
5    Square pyramidal, trigonal bipyramidal, pentagonal planar
6    Octahedral, trigonal prismatic, pentagonal pyramidal
7    Pentagonal bipyramidal
8    Square antiprismatic
==== ====

For every symmetry, the following information is stored in a manner that permits
compile-time computation:

- Ligand to ligand angle function
- A minimal set of rotation mappings that allow the generation of all
  superimposable position sequences
- A set of tetrahedron definitions whose signed volumes permit the distinction
  of chiral differences between position sequences


Stereopermutations
==================
TODO Text about the symbolic calculation of stereopermutations
In many cases, the relative positioning of a central atom's substituents
corresponds to a common local symmetry. Most of organic chemistry, for instance,
can be composed of linear, bent, trigonal planar and tetrahedral local
symmetries. The most common source of chiral information in organic chemistry is
an asymmetric tetrahedral carbon, where its four subtituents are all different
(i.e. have different ranking according to the IUPAC rules). From a local
symmetry and the information that all substituents are distinct
ranking-wise, we can conclude that there are two non-superimposable
permutations of ligands. Stereopermutations of any symmetry can be computed
symbolically[REF], including cases in which substituents are mutually linked,
such as in multidentate ligands. Special care must be taken in order to reduce
haptically bonded ligands to the correct local symmetries and ranking. In
addition, symbolic computation of permutations does not consider the conditional
feasibility of assembling linked ligands depending on a symmetry's
ligand-to-ligand angle, bridge lengths, element types and bond types involved in
the link.


Molecule model
==============

The molecule data model consists of two parts: A graph, and a list of
stereocenters.


Graph
-----
Molecules are modeled as an undirected graph consisting of atom vertices and
bond edges. Vertices store the atomic element type while vertices store a bond
type that distinguishes the bond orders one through six as well as an aromatic
and a so-called eta bond, which models connections between a central atom and a
haptically-bonded subset of atoms (i.e. a contiguous group of atoms within a
separate molecule that contribute to the formation of the haptic bond).

A molecule's graph must consist of a single connected component, meaning that
there must be a path from any atom of the molecule to any other. Single atoms
are not considered molecules, so a molecule must consist of at least two
mutually bonded atoms. Removing atoms or bonds from diatomic molecules are
disallowed operations. Similarly, disconnecting a molecule into two logical
molecules by removing a particular bond or atom is also disallowed.

Every non-terminal atom (i.e. any vertex with more than one edge) may carry
geometric and stereocenter information. 


Stereocenters
-------------
Currently, this library models two types of stereocenters: The first models the
atom-centric relative placement of substituents across geometries ranging from
linear to square-antiprismatic. The second type models differences owing to
strong bond-centric rotational barriers at double bonds.

The two stereocenter types cannot overlap. While an atom-centric stereocenter
models an idealized symmetry at a singular atom, a bond-centric stereocenter
models the individual geometries at both endpoint atoms.

In order to determine whether a particular atom in a molecule is an atom-centric
stereocenter, its substituents are ranked according to an implementation of the
IUPAC sequence rules as laid out in the 2013 Blue Book. The rankings are
transformed into an abstract substituent case (e.g. octahedral (A-A)BBCD) and a
symbolic computation is carried out to determine the number of permutations that
are not superimposable via spatial rotations within the idealized symmetry. The
set of resulting permutations of the substituent symbols is called the set of
stereopermutations. If this set contains more than one stereopermutations, then
the atom is an atom-centric stereocenter under that symmetry.

If substituents are haptic or multidentate, an additional algorithm removes
stereopermutations it deems obviously impossible. All bridge lengths between
pairs of chelating atoms of a multidentate ligand are checked against the atom
pair's bite angle within the idealized symmetry. Additionally, haptic ligands'
cones are checked to ensure they do not overlap. Indices within the set of not
obviously impossible stereopermutations are called assignments.

Bond-centric stereocenters are less complex, since only simple cis/trans
isomerism is currently modeled. If a prospective bond is a double bond and its
endpoint atoms both have three substituents, then both atoms' substituents are
ranked. If neither endpoints' two outward-facing substituents rank equally, then
the bond and its two endpoints are considered a bond-centric stereocenter. This
model of rotational energy barriers is deficient in a multitude of cases (see
Known Issues).
