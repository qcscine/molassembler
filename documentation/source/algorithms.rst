=====================
Underlying Algorithms
=====================

.. role:: strikethrough
   :class: strikethrough

.. role:: green
   :class: green

The aims and behavior of several algorithms must be clearly understood to avoid
pitfalls and misconceptions about their applicability.


Bond discretization
===================
The discretization of fractional bond orders to bond types is currently a fairly
boorish round-to-nearest-integer algorithm and is hence particularly error prone
for particularly weakly-bound metal ligands (around 0.5) and aromatic bonds
(around 1.5). Until conjugated system detection is implemented and this
algorithm is improved, you may want to consider using binary bond discretization
(bonded / not bonded) if you work with conjugated systems.


Symmetry determination
======================
If no coordinates are present when constructing a Molecule, the idealized
symmetry on any non-terminal atom must be determined from its local
substituents within the set of appropriately sized symmetries. For main-group
atoms, this is accomplished by application of a very basic valence shell
electron pair repulsion (VSEPR) algorithm. No symmetry determination algorithms
are currently implemented for non-main-group atoms, and the first symmetry of
appropriate size is chosen.


Ranking
=======
In order to establish the relative priority of an atom's substituents, a ranking
algorithm is applied that follows the IUPAC recommendations laid out in the 2013
Blue Book [1]_, generalized to larger symmetries. The IUPAC sequence rules are as
follows, slightly reworded:

=== =============
Nr. Sequence rule
=== =============
1a  Higher atomic number precedes lower. Atomic numbers in Mancude rings and
    ring systems are the mean of the atomic numbers in all mesomeric Kekule
    structures.
--- -------------
1b  A duplicate atom node whose corresponding nonduplicated atom node is closer
    to the root ranks higher than a duplicate atom node whose corresponding
    nonduplicated atom node is farther from the root
--- -------------
2   Higher atomic mass number precedes lower
--- -------------
3   When considering double bonds and planar tetraligand atoms Z precedes E and
    this precedes nonstereogenic double bonds.
--- -------------
4a  Chiral stereogenic units precede pseudoasymmetric stereogenic units and
    these precede nonstereogenic units
--- -------------
4b  When two ligands have different descriptor pairs, then the one with the
    first chosen *like* descriptor pairs has priority over the one with a
    corresponding *unlike* descriptor pair. Descriptors are alike if they are
    within the same set. The two sets are {R, M, Z} and {S, P, E}.
--- -------------
4c  r precedes s and m precedes p
--- -------------
5   An atom or group with descriptors {R, M, Z} has priority over its
    enantiomorph {S, P, E}
=== =============

The following alterations arise due to the current implementation state:

- Stereodescriptors are transformed into the index of permutation within the
  symmetry's abstract ligand case set of stereopermutations. The descriptors R/S
  correspond to indices of permutation within the set of stereopermutations of
  the abstract tetrahedral ABCD ligand case. Z/E exist as indices of
  permutations within the set of stereopermutations of the bond-centric
  stereocenter. M and P stereocenters are not modelled.
- For larger or smaller symmetries, sequence rule 5 is altered to give priority
  according to the index of permutation: {R = 1, Z = 1} > {S = 0, E = 0}.
- *Like* is altered to mean stereodescriptors with the same index of
  permutation.
- No distinction is made whether a stereocenter is merely psuedoasymmetric. When
  based on indices of permutation, sequence rules 4c and 5 can be conflated, and
  I cannot think of a molecule in which removing the distinction changes the
  sequence rule application for sequence rule 4a. Additionally, 
  it does not affect whether the exact chirality of that stereocenter must be
  fixed in the final conformation for this library. Keeping this information
  merely permits direct generation of the enantiomer for molecules in which
  each stereocenter has only two stereopermutations. Since we aim to be
  applicable to stereocenters with many more stereopermutations, the value of
  the distinction is slight. Sequence rule 4a is shortened to enforce priority
  of stereogenic units over nonstereogenic units. Sequence rule 4c is conflated
  with 5.
- No functionality currently exists to alter atomic mass number. Sequence rule 2
  is discarded.
- Conjugation detection is currently not implemented. No regularization of
  Mancude rings and rings systems is performed. This is the only change that can
  lead to false differences between substituents when they are actually
  identical.

In summary, the applied sequence rules are:

=== =============
Nr. Sequence rule
=== =============
1a  Higher atomic number precedes lower. :strikethrough:`Atomic numbers in
    Mancude rings and ring systems are the mean of the atomic numbers in all
    mesomeric Kekule structures.`
--- -------------
1b  A duplicate atom node whose corresponding nonduplicated atom node is closer
    to the root ranks higher than a duplicate atom node whose corresponding
    nonduplicated atom node is farther from the root
--- -------------
2   :strikethrough:`Higher atomic mass number precedes lower`
--- -------------
3   When considering double bonds and planar tetraligand atoms Z precedes E and
    this precedes nonstereogenic double bonds.
--- -------------
4a  Chiral stereogenic units precede :strikethrough:`pseudoasymmetric
    stereogenic units and these precede` nonstereogenic units
--- -------------
4b  When two ligands have different descriptor pairs, then the one with the
    first chosen *like* descriptor pairs has priority over the one with a
    corresponding *unlike* descriptor pair. Descriptors are alike if they
    :green:`have an equal number of stereopermutations and equal index of
    permutation` :strikethrough:`are within the same set. The two sets are {R,
    M, Z} and {S, P, E}`.
--- -------------
4c  :strikethrough:`r precedes s and m precedes p`
--- -------------
5   An atom or group with :green:`higher index of permutation has priority`
    :strikethrough:`descriptors {R, M, Z} has priority over its enantiomorph {S,
    P, E}` 
=== =============


Cycle detection
===============
Cycle detection and enumeration is handled by the excellent library
RingDecomposerLib [2]_, which avoids exponential space requirements using Unique
Ring Families [3]_.

Conformer generation
====================
For each sought conformation, unassigned stereocenters in the input molecule are
assigned randomly according to the relative statistical occurrence weights of
the stereopermutations.

Conformer generation itself is based on four-dimensional Distance Geometry [4]_.
This library's implementation features the following:

1. A spatial model generates atom-pairwise bounds on their distance in the final
   conformations and four-atom chiral constraints when distance bounds cannot
   encompass chiral elements of complex symmetries. For large symmetries, chiral
   information is captured by using multiple chiral constraints.
2. The distance bounds are smoothed to conform to the triangle inequalities.
   After each successive choice of a fixed distance between the bounds, you can
   choose to re-smooth all bounds (full metrization) or stop re-smoothing after
   a fixed number of chosen distances (partial metrization). Smoothing is
   performed by transferring the problem to a graph shortest-paths problem [5]_
   and finding the shortest paths with the GOR1 algorithm [6]_ instead of the
   naive Floyd-Warshall algorithm.
3. The bounds are embedded in four dimensions and refined in two stages,
   permitting the chiral constraints to invert by expanding into four
   dimensions, and then compressing the fourth dimension back out.
4. The refinement error function is modified to enable the placement of haptic
   ligand's bonding atoms' average position at symmetries' idealized ligand
   sites.


References
==========
.. [1] Favre, H.A., Powell W.H. Nomenclature of Organic Chemistry: IUPAC
   recommendations and preferred names. Royal Society of Chemistry. **2013**.

.. [2] Flachsenberg, F.; Andresen, N.; Rarey, M. RingDecomposerLib: An
   Open-Source implementation of Unique Ring Families and Other Cycle Bases. *J.
   Chem. Inf.  Model.*, **2017**, 57 (2), pp 122–126

.. [3] Kolodzik, A.; Urbaczek, S.; Rarey, M. Unique Ring Families: A Chemically
   Meaningful Description of Molecular Ring Topologies. *J. Chem. Inf. Model.*,
   **2012**, 52 (8), pp 2013–2021

.. [4] Blaney, J.M.; Dixon, J.S. Distance Geometry in Molecular Modeling. *Rev.
   Comp. Ch.* **2007**, pp. 299-355

.. [5] Havel, T.; Wüthrich. K. A distance geometry program for determining the
   structures of small proteins and other macromolecules from nuclear magnetic
   resonance measuremenets of intramolecular 1H-1H proximities in solution *B.
   Math. Biol.* **1984**, 46.4, 673-698.

.. [6] Cherkassky, B. V., Goldberg, A. V., & Radzik, T. Shortest paths
   algorithms: Theory and experimental evaluation. *Math. Prog.*, **1996**.
   73(2), 129–174.
