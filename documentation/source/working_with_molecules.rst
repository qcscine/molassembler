======================
Working with molecules
======================


Molecule construction
---------------------

.. graphviz:: graphs/information_flow.dot

In order to create a molecule instance, information can be supplied and inferred
in many ways. 

Minimally, the atoms' element types and their connectivity are needed to
constitute a molecule's graph. From a graph, the presence of stereocenters can
be inferred by using algorithms to determine the local symmetry at an atom using
information present in the graph about its substituents. Any inferred
setereocenters are unassigned, i.e. generated conformers will be composed of
all possible stereopermutations at those stereocenters.

Atom connectivity can also be discretized from a bond order matrix, which in
turn can be very roughly approximated using pairwise atom distances from spatial
information.

Atom-centric symmetries and stereocenter assignments can be gleaned from spatial
information.

Molecular file formats supply varying levels of information about molecules,
ranging from the XYZ format, which supplies merely element types and spatial
information, to the MOLFile format, which supplies element types, discretized
connectivity and spatial information.


Options
-------
The global behavior of the library can be altered in several respects.

Whether an atom is a stereocenter is, in some cases, dependent on ambient
temperature. Trigonal-pyramidal nitrogen can invert rapidly at room temperature,
but is a stereocenter at low temperature. By default, the library is in the
high-temperature approximation, where trigonal-pyramidal nitrogen centers are
can only be a stereocenter if it is part of a small cycle (sizes 3, 4).

While editing a molecule on the graph level, the library attempts to propagate
existing chiral state in its list of stereocenters to the new graph. By default,
state is propagated only if a singular mapping of indices to the new symmetry
exists in which no symmetry angles are altered. Which propagations of chiral
state are acceptable can be altered in the header `Options.h`.


Editing
-------
The maxim of the library is to permit the specification of any desired graph,
however unreasonable, within the logical constraints of the molecular model.
These constraints are that a Molecule class can only contain a single connected
component and must consist of at least two atoms. Apart from that, you can
specify really weird graphs if you so desire. Keep in mind, however, that
particularly weird graphs can be at odds with either the reduction of
stereopermutations to assignments, where algorithms throw out obviously
impossible permutations, or with the library's spatial modeling.


Conformer generation
--------------------

- Symmetry idealization strictness
