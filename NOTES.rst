CONTINUE AT
-----------
- Metrization during distance matrix generation in DistanceBoundsMatrix
  (At step 7 of DG steps from p.15)

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
- 4D DGRefinement: gradient does not yet encompass compress=true, gradient +=
  2*(4th dimension)
- Sanity test should encompass EVERY possible assignment being generated AND
  recovered! Remember that currently, the alternate set of tetrahedra
  definitions is in use!
  - Implementation is written, but IterateStereocenterPermutations is in dire
    need of some tests before anything can be truly verified
- Some tests are still more like analysis scripts, and some analysis scripts
  could be reformulated as tests
- There are now four different implementations of the DG process: In
  generateConformation, and spread across tests/ and analysis/. As soon as the 
  main one is reliable, refactor the rest to the main one.
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
