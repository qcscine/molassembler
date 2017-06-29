TODO
----
- CNStereocenter and EZStereocenter will clash! Grubbs cat::

              R R     R
               \|    /
    sq. py. ->  W = C
               /|    \
              R R     H

  CNStereocenter on wolfram, an EZStereocenter on W and C. How to handle this?

  Is it acceptable to just instantiate the CNStereocenter on W in this case?
  There is missing dihedral information, but as long as C is recognized as
  trigonal planar the overall structure will probably be correct! The twist
  around the W=C bond is probably restricted, but is it clear which conformation
  is best?

- Metrization during distance matrix generation in DistanceBoundsMatrix
  (At step 7 of DG steps from p.15)
- Parallelize DG
- Make dependence on alternate set of Symmetry tetraeder subdivisions clearer
  in code / analysis
- Set up CMake properly
- Consider penalization of unsuitable geometries (using VSEPR /
  determineLocalGeometry to judge) in SymmetryFit -> see what happens when you
  try to load testosterone, which has some strained tetrahedral centers that are
  fitted as seesaws!
- Investigate link-time optimization
- Make sure EZStereocenter is stable against the situation where (and the twist
  is a given) -> Also, the involvedAtoms of this case overlap! No singular
  stereocenter per atom in the entire molecule! Unless you create another type
  that handles this case specifically, involving all three atoms.::
    
    1
     \
      3 = 4 = 5 - (67)
     /
    2

- Rewrite the ranking algorithm as a BFSVisitor
- The Readers really need a file exists check
- Fixed atoms in DG -> how? Probably by just specifying the geometry without
  variance in the distance bounds and then setting the fixed atoms' gradients to
  zero in the refinement stage
- Transition to property-based testing with rapidcheck (github) and/or fuzz
  testing
- Atom removal safety of code -> getNumAtoms, getNumBonds, etc. Make the full
  set of data be contiguous every time (atom indices range from 0 -> nAtoms - 1
  AdjacencyList may also be prone to errors in this regard
- Aromatic cycles have to be detected for DG! This property has to be detected
  for generated structures to be more reasonable, but also not necessarily from
  the start. But watch out, don't forget to integrate it at some point! A
  starting point is currently in include/repurpose/
- Use LFT to determine which geometry? MO-level calculations, perhaps
  approximable with low cost
- Give some thought to how to treat charge / total electron counts
- Make everything nothrow as best as possible
- Reconsider float <-> double necessity, especially in DG: Might get some
  significant gains by using floats
