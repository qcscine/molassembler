CONTINUE AT
-----------
- How to get angles between atoms in graph? Must include some conception of
  charge
- Constraint collection and population of DistanceBoundsMatrix
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

4. Implement a function that permutes all CN4Stereocenter instances
5. Implement basic DG
6. Demonstrate functionality with a very simple example, e.g. CH(Cl)(Br)(I)
   MOLFile, then permute and generate 3D structures of both stereoisomers.

- Fixed atoms in DG -> how?
- Transition to property-based testing with rapidcheck (github) and/or fuzz
  testing
- Atom removal safety of code -> getNumAtoms, getNumBonds, etc. Make the full
  set of data be contiguous every time (atom indices range from 0 -> nAtoms - 1
  AdjacencyList may also be prone to errors in this regard
- Add hooks to git to automatically build a release version on commit and run
  the tests
- Should PositionCollection really be a member of Molecule? I don't think so
  Perhaps optionally, or better yet, cached
- Should AromaticRing really be a GraphFeature? Isn't that somewhat a misnomer
  anyway? The whole necessity for their existence was that the connectivity of
  vertices and edges is sometimes insufficient to fully specify a molecule's
  shape. Is an aromatic ring such a thing? The property of specific local
  connectivities on specific types of atoms + a cycle directly leads to the
  planarity of the involved atoms. It's not something that HAS to be specified
  separately for those atoms. Separately, for DG, this special property has to
  be detected for generated structures to be more reasonable, but also not
  necessarily from the start. But watch out, don't forget to integrate it at
  some point! It's currently in include/repurpose/
- IO.h: will have to be changed eventually to call DG to generate a 3D
  structure if there is none.  Maybe cache 3D structures? Additionally,
  modernize it to use C++17's filesystem TS
- Use LFT to determine which geometry? MO-level calculations, perhaps
  approximable with low cost


Is atomic charge fully determined?
----------------------------------

The information present is the atom type of the center atom, the atom types of
all neighbors, and the bond types.


Interface
---------

- How do I want to conceptualize a Molecule, and how should working with one
  look in practice?

  Current idea::

  // Constructors
  // 1 reads file, creates adjacency list if file type does not contain
  // connectivity information
  Molecule init_from_file(Filename("asdf.xyz")); 
  // 2 from Delib AtomSet, detects adjacency list
  Delib::AtomSet atom_set();
  Molecule init_from_AtomSet(atom_set); 
  // 3 copy constructor
  Molecule copy_from_other(init);

  // Permutators
  vector<Molecule> structures = {
      init_from_file, 
      init_from_AtomSet,
      copy_from_other
  };

  structures = permute_hydrogen_replacements(structures);
  
  // code in permutator
  vector<Molecule> permute_*(const vector<Molecule>& structures) {
      vector<Molecule> permutated;
      for(const auto& molecule : structures) {
      }
  }

