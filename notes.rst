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
