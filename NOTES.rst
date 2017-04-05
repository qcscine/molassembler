CONTINUE AT
-----------
- Metrization during distance matrix generation in DistanceBoundsMatrix
  (At step 7 of DG steps from p.15)
- Minimization in generateConfiguration.h

Things that need tests
----------------------

- AdjacencyListAlgorithms
- BondDistance
- CN4Stereocenter
- Cache
- CommonTrig
- IO
- Molecule (when finished)
- StdlibTypeAlgorithms
- StereocenterList
- TreeAlgorithms
- DistanceGeometry
  
  - DistanceBoundsMatrix
  - MetricMatrix
  - generateConformation (when finished)


TODO
----
- The recent refactor of the symmetry fitting and separation into analysis and
  testing binaries broke the chirality-constraints testing due to a change of
  format. The entire application needs a good logger that can selectively output
  specific information. Then you can go about fixing the symmetry fit output
  from the adjacencylist and then fix the regression in the data format for the
  R analysis. Is it really best to add a logger to get specific output at
  selected points in the program?
  
  Logger:
  + No changes to overall structure
  - Inefficiency due to logging level checks

  Change structure: (Get Positions, AdjacencyList from IO, don't immediately
  construct a molecule, pass a parameter to AdjacencyList's inferFromPositions
  method)
  - Small change structurally, adds a singular-purpose IO interface
  + No need for a logger, no added inefficiency
  

- Make sure EZStereocenter is stable against the situation where (and the twist
  is a given) -> Also, the involvedAtoms of this case overlap! No singular
  stereocenter per atom in the entire molecule! Unless you create another type
  that handles this case specifically, involving all three atoms.
    
    1
     \
      3 = 4 = 5 - (67)
     /
    2

- To create a 3D structure that needs to create a specific assignment on a
  stereocenter, it might be preferable to set incorrect and unique 1-2 distances
  to the central atom during coordinate generation (with 1-3 distances, this
  enforces the creation of the correct stereocenter) and then pay the price
  during minimization (where we correct the 1-2 distances) instead of generating
  the wrong stereocenter with correct 1-2 distances and then going through
  high-energy contortions to correct stereocenters to the right chirality.
  Another idea would be to contract all falsely expanded bonds manually
  immediately after generating the correct stereocenter, but the issue with this
  are cycles. Any modifications in bond lengths would negatively affect cycle
  conformation sampling. Although I would expect this is the same when a
  stereocenter that involves a cycle has to be corrected It would be straight up
  amazing though - we could potentially avoid chirality constraints altogether
  during conformational generation with the issue that we need them anyway when
  we try to detect stereocenters (although there we could mess with the
  structures too, forcing 1-3 deviations). But it doesn't allow me to just
  straight up ignore chirality constraints. I have to do both implementations so
  I have something to compare.
- Maybe some functionality of the EZStereocenter / CNStereocenter items can be
  moved to the symmetry library. It would be a handy place to keep all the
  important algorithms of extracting constraints, fitting to symmetries /
  assignments, etc. But maybe it's too tightly coupled with AdjacencyList...
- Rewrite the ranking algorithm as a BFSVisitor
- 4D DG is still missing
- The Readers really need a file exists check
- Fixed atoms in DG -> how? Probably by just specifying the geometry without
  variance in the distance bounds and then setting the fixed atoms' gradients to
  zero in the refinement stage
- Transition to property-based testing with rapidcheck (github) and/or fuzz
  testing
- Atom removal safety of code -> getNumAtoms, getNumBonds, etc. Make the full
  set of data be contiguous every time (atom indices range from 0 -> nAtoms - 1
  AdjacencyList may also be prone to errors in this regard
- Add hooks to git to automatically build a release version on commit and run
  the tests
- Aromatic cycles have to be detected for DG! This property has to be detected
  for generated structures to be more reasonable, but also not necessarily from
  the start. But watch out, don't forget to integrate it at some point! A
  starting point is currently in include/repurpose/
- Use LFT to determine which geometry? MO-level calculations, perhaps
  approximable with low cost


Is atomic charge fully determined?
----------------------------------

The information present is the atom type of the center atom, the atom types of
all neighbors, and the bond types.
