# In-progress notes
- More granular CNStereocenter and PermutationState tests needed
- Cone overlap check missing
- Eta bond dynamism is missing
- DG visualization needs a replacement for tetrahedronHighlights

# Haptic ligands implementation plan

## Approaches
1. Generalize CNStereocenter to deal with haptic ligands.
   
   + No new class, possibly reduces graph handling complexity
   - Stereocenters without haptic ligands (wide majority) have additional
     overhead due to additional abstraction

2. Add a new Stereocenter class that can handle haptic ligands, internally
   wrapping the functionality of CNStereocenter.

## Eta bond type
- For every metal atom in the molecule, check if any adjacents are pairwise
  adjacent to each other. Groups of contiguous adjacents are to be considered
  as a ligand to the metal, and their internal bond type to the metal marked
  Eta.
- Stereocenters on non-metal ligand atoms are to disregard eta bonds in
  determining which atoms are adjacent, fitting geometries, et cetera.
- Stereocenters on metal ligand atoms group eta bonds as singular ligands
- Stereocenters on metal ligand atoms must consider eta bonds in all regards.

## Class handling haptic bonds specification
- The final class accepts, as construction arguments,
  - the regular ranking result of all adjacent atoms
  - a nested vector detailing which indices are grouped as a ligand
  - Links between its atom indices in the graph
- Sanity checks to ensure that impossible assignments are removed before their
  three-dimensional representation is attempted:
  - Any assignments in which the closest distance between unlinked haptic
    ligand cones is smaller than the sum of the respective smallest atom vdw
    radii -> Risk of cone intersection
- Coordinate fitting to decide assignment and geometry: Use haptic atom centroid
  as cone base center.

## DG implications
- Modeling of haptic ligands:
  - The subset of haptically bonded atoms of the overall ligand are typically
    coplanar. In a first approximation, we model them as cyclic polygons.
  - Together with the metal center, haptic ligand atoms are placed on the base
    circle of a cone. The base radius is the circumradius of the cyclic polygon
    constructed using the bond lengths of the haptic ligand bonds. The height
    of the cone is an average of the individual element type's eta bond
    distance to the metal type.
  - If the haptic ligand consists of two atoms only, the base radius is exactly
    half the bond distance between them.
  - Angle bounds between individual atoms on two ligands, which may possibly
    both be haptic, are calculated as the underlying symmetry's ligand angle
    plus or minus both ligand cones' apex angles.
  
- Usage of angle(i, j, k) where i, j, k are atom indices is no longer viable in
  DG together with an absolute variance parameter for those stereocenters that
  carry a haptic ligand.

## Inadequacies
- It is probably slightly inaccurate to model haptic ligands as right circular
  cones whose base center is the circumcenter of the cyclic polygon. In the case
  of significant bond length variance in the polygon edges, the circumcenter may
  even lie outside the polygon. Particularly in those cases, it is probably more
  accurate to estimate that the metal atom is placed perpendicular to the
  atomic centroid. It is not clear how to adapt the cone intersection test and
  individual angle calculations to this possibly improved approximation,
  however. Additionally, considering that it is considerably rare in the
  literature to encounter polyhapto heterocyclic ligands, it is unlikely that
  such distorted polygons will be encountered where the approximation is
  problematic.
