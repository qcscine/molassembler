one issue in DG that can probably be improved upon is that for cycles smaller
than 6, meaning 3, 4 and 5-sized cycles, the overall geometries for the
individual atoms in the cycle are distorted away from the ones determined by the
graph.

The edge lengths are known, but that is all the information we get for all
polygons. 

- How should angles to ring atom substituents be calculated, depending on how
  many atoms are there? what happens if e.g.::
      
         .·C
    R ~ W  |
         °·C

  where wolfram could have *any* coordination geometry?

- Triangles are always flat, so for cycles of size 3 we can calculate all
  internal angles using the law of cosines.
- There are no aromatic cycles of size 4 (I think), so there is no particular
  use in calculating the angles properly and adding variance to the purportedly
  ideal geometry
- There are, however, aromatic cycles of size 5, and it is probably ideal to
  maximize their area. See if you can find examples of converged aromatic
  5-cycles, including heavily hetero-substituted ones, and see if they are
  concyclic!
- There exist formulas to determine the maximum pentagon area from the edge
  lengths alone (recent papers), and also the circumradius. From that, it should
  be easy to construct a pentagon and get the internal angles.
- Maybe it is still better to just ignore this (admittedly minor) optimization
  and just increase tolerance to the strained cycles to account for distortion.
  However, I'm pretty sure the deviation from the ideal geometry in strained
  cycle is does not have the same angle averages, so it would be more accurate
  to properly model the cycles using geometric considerations.
- Maybe it is best to investigate the internal angle distributions for 3, 4 and
  5-cycles to ensure the best path going forwards, perhaps using converged
  calculations of molecules from the CCCDB... But no way to ensure these
  considerations hold for any heavily substituted cycles!
