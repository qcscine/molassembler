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

It is recommended to fix the connectivity of a molecule with previously known
connectivity information. The discretization of bond orders to connectivity
currently has pitfalls for conjugated systems, and the calculation of bond
orders from spatial distances is fragile.

The best supported file format for the extraction of assignments of
stereocenters is therefore the MOLFile.


Options
-------
The global behavior of the library can be altered in several respects.

In some cases, whether an atom is a stereocenter is dependent on ambient
temperature. Trigonal-pyramidal nitrogen can invert rapidly at room temperature,
but is a stereocenter at low temperature. By default, the library invokes the
high-temperature approximation, where trigonal-pyramidal nitrogen centers can
only be a stereocenter if it is part of a small cycle (sizes 3, 4) where the
strain hinders inversion.

While editing a molecule on the graph level, the library attempts to propagate
existing chiral state in its list of stereocenters to the new graph. By default,
state is propagated only if a singular mapping of indices to the new symmetry
exists in which no symmetry angles are altered. Which propagations of chiral
state are acceptable can be altered.

The pseudo-random number generator is centralized and can be re-seeded for
reproducible behavior.


Editing
-------

The maxim of the library is to permit the specification of any desired graph,
however unreasonable, within the logical constraints of the molecular model.
These constraints are that a Molecule class can only contain a single connected
component and must consist of at least two atoms. Apart from that, you can
create completely unreasonable graphs. Keep in mind, however, that these can be
at odds with either the reduction of stereopermutations to assignments, where
algorithms throw out obviously impossible permutations (yielding stereocenters
with no possible assignments) or with the library's spatial modeling in
conformer generation.


Conformer generation
--------------------

The generation of conformational ensembles of molecular graphs has a few
peculiarities you need to know about.

For one, it is very explicit in the library that conformer generation can fail.
The possibility of failure is expressed by the explicit Result type of conformer
generation interfaces, which is essentially a two-state variant type holding
either your requested conformational ensemble or some information about the
reason for failure. The possibility of failure is explicit precisely to force
you to consider it when integrating this library. A molecule you have generated
may have stereocenters with no possible assignments due to too short bridge
lengths for multidentate ligands or too many large haptic ligands whose cones
overlap in any arrangements. It may also conflict with the spatial model, which
has upper limits on granted slack for particularly strained graphs.

Furthermore, the conformer generation sequence will increase modeling tolerances
for each progressive embedding failure across the entire molecule. Even if only
one part of your molecule is strained and requires this additional effort, the
spatial model for the entire molecule will slacken, reducing local symmetry
idealization strictness even in clear-cut cases. Heavily loosened tolerances can
lead to odd-looking spatial structures. In those cases, you may want to further
refine these with force field or semiempirical quantum chemical methods before
running any more precise calculations.
