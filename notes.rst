Files
-----

- I_MoleculeCollection.hpp

  An abstract base class defining the basic interface for a MolecularGraph
  class. This should store pointers to Molecule instances and hand them out to
  classes that want to work with them.

- Molecule.hpp

  A struct that will contain all info required of a molecule. Only few
  algorithms, operators will be contained here

- MoleculeManip.hpp

  Examples of functions that can manipulate collections of molecules to
  structure generation ends and a sample implementation of folding them onto a
  set of initial structures.

  Thought as a start to re-implement the structure generation project.


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



