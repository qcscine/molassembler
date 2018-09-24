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
stereocenters are unassigned, i.e. generated conformers will be composed of
all possible stereopermutations at those stereocenters.

Atom connectivity can also be discretized from a bond order matrix, which in
turn can be very roughly approximated using pairwise atom distances from spatial
information.

Atom-centric symmetries and its stereocenter assignments can be gleaned from
spatial information.

Bond-centric stereocenters are interpreted at each bond where both ends have an
assigned atom-centric stereocenter. If the supplied conformation matches a
permutation of the bond-centric stereocenter, the bond-centric stereocenter is
assigned and kept. Incidental matches and undesired dihedral freezes can be
avoided by specifying a minimal bond order for which a rotational barrier can
exist. Any bonds whose bond orders do not meet this threshold are not considered
for bond-centric stereocenter interpretation.

Molecular file formats supply varying levels of information about molecules,
ranging from the XYZ format, which supplies merely element types and spatial
information, to the MOLFile format, which supplies element types, discretized
connectivity and spatial information.

Recommendations
---------------
The preferred method for mass molecule construction is binary discretization
from supplied fractional bond orders and coordinate information.

Nearest-integer bond discretization has pitfalls for conjugated systems, and the
calculation of bond orders from spatial distances only is unreliable due to the
limited underlying model of UFF bond distances.


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

Information about a Molecule can be obtained from its constituting graph and
stereocenter list classes, but editing is centralized at the Molecule interface
since graph and stereocenter changes force a full re-ranking of all atoms and
possibly changes in the list of stereocenters.

Among the possible changes in the list of stereocenters on re-ranking are:

- Stereocenters may change the number of possible assignments. This means a
  stereocenter may become chiral (i.e. there are multiple stereopermutations) or
  cease to be chiral (there is only one stereopermutation).
- If a stereocenter is assigned, its assignment may change.

It is therefore important that **upon each molecular edit**, any state about or
within the StereocenterList of the molecule being edit that is stored externally
to the instance itself is considered invalidated.


The Eta bond type
-----------------

On a graph level, each constituting atom of a haptic ligand that is part of the
bond of the ligand to a transition metal must have a bond to that transition
metal. These bonds are of the Eta bond type, regardless of their individual
order. It is, however, the task of the library to differentiate which bonds
between atoms are to be considered an Eta bond type and which are not, based
wholly upon the graph structure. The public API of molassembler does not accept
user-defined Eta bond types. Specify Single as bond type in the creation of
haptic bonding.


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
