Remaining DG deficiencies
-------------------------
- Strain imposed on EZStereocenter substituents that are in small cycles is not
  accounted for::
            
          . C
    O = C   |
          ° C

  EZStereocenter's angle(C, C, C) returns 120°, and the resulting bounds are not
  modified if the involved atoms are cycle members
- Cycles of size 8 can have considerable strain due to triple bonds in the
  cycle, this is not considered
- Bounds are perhaps a little too loose now, spiro centers look awful
- Double bonds in DG maybe shouldn't enforce absolute flatness, but have
  tolerance. Either that or reduce the error function contribution of
  Flat-target chirality constraints by a factor to leave 1-2 distance bounds
  unviolated in strained molecules (see strained-db-aromatic-multicycle
  examples)
           
TODO
----
- Recheck geometry_assignments relative weights for correctness
- Fix all the failing tests
- Stereocenters that merely mark non-standard geometries and do not provide any
  assignments are probably not included in ranking at all, but they should be!
- Documentation:
  - README.md file
- EZStereocenters in small cycles aren't truly stereocenters. There's no reason
  to be permuting 3 EZStereocenters in a benzene, any E arrangements are
  infulfillable anyway. But in oct-1-ene, E/Z differences are realizable. How
  should this be approached?
  IterateStereocenterPermutations should account for that.
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

  This is the approach taken now, but no tests to ensure this is handled
  correctly.

- Make dependence on alternate set of Symmetry tetraeder subdivisions clearer
  in code / analysis
- Make sure EZStereocenter is stable against the situation where (and the twist
  is a given) -> Also, the involvedAtoms of this case overlap! No singular
  stereocenter per atom in the entire molecule! Unless you create another type
  that handles this case specifically, involving all three atoms.::
    
    1
     \
      3 = 4 = 5 - (67)
     /
    2

- The Readers really need a file exists check
- Fixed atoms in DG -> how? Probably by just specifying the geometry without
  variance in the distance bounds and then setting the fixed atoms' gradients to
  zero in the refinement stage
- Aromatic cycles have to be detected for DG! This property has to be detected
  for generated structures to be more reasonable, but also not necessarily from
  the start. But watch out, don't forget to integrate it at some point! A
  starting point is currently in include/repurpose/
- Use LFT to determine which geometry? MO-level calculations, perhaps
  approximable with low cost
- Give some thought to how to treat charge / total electron counts
- Since DG is stochastic, conformational ensemble generation is not systematic.
  It would be useful to consider if there are simple ways of manipulating DG to
  get systematic exploration of rotatable bonds. An alternative to that is at
  least a measure that quantifies how much of the conformational space has been
  sampled with a certain number of generated conformers. There are two ways of
  going about this probably, one being a statistical consideration of the
  intervals of allowed dihedrals, and the other being an estimation of the
  number of conformers based on the two symmetries at opposite ends of the
  rotatable bond and then multiplying all these together. It is probably
  impossible to consider that choosing one dihedral may limit other dihedral
  angles due to vdw radius overlap, and consequentially it is probably necessary
  to consider all dihedral angles as independent. Single bonds in small cycles
  do not have the same range of possible values as ones in linear sections, and
  the same applies to aromatic systems. This has to be considered in rotamer
  counting.

Improvement considerations
--------------------------
- NL.16 Use a conventional class member declaration order
  types -> constructors, assignments, destructor -> functions -> data
  (public - protected - private)
- Consider replacing std::set with set::unordered_set where no in-order
  traversal or lexicographical comparison needed
- Replace boost::optional<bool> with boost::tribool
- Set up CMake properly
- Investigate link-time optimization
- Parallelize DG (?)
- Transition to property-based testing with rapidcheck (github) and/or fuzz
  testing
- Make everything nothrow as best as possible
- Reconsider float <-> double necessity, especially in DG: Might get some
  significant gains by using floats
- Consider a class that provides the same interface as Molecule, but, instead
  of losing chiral information on e.g. a ligand addition when the chiral
  information preservation mode is on EffortlessAndUnique and we're going from
  square planar to square pyramidal, it just internally splits into the two
  possible assignment cases for that molecule. Modifications are only valid if
  applying them on both Molecules doesn't throw. That would be a cool way to
  keep complete chiral state during modifications. This may be implementable
  with an additional CNStereocenter function that returns the possible
  assignments in the target symmetry without altering internal state.
  Then you can alter with a lower preservation mode, copy and assign as desired.
  This would allow free-standing functions with signatures like::

    std::vector<Molecule> addBond(const Molecule& mol);

General notes
-------------
- License for URF library?
- Various strained organic molecules are taken from "Survey of strained org
  molecules" by Liebman, Greenberg. 1975
