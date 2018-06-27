============
Known Issues
============

Ranking algorithm deficiencies
------------------------------
No regularization of Mancude rings and rings systems is performed in sequence
rule 1. This can lead to false differences between substituents when they are
actually identical.

Rotational isomery
------------------

The current data model of bond-centric stereocenters can treat a large number of
cases correctly, but is insufficient. For instance, allene systems' rotational
isomery cannot be captured with the current model. The current limitation to
trigonal pyramidal - double bond - trigonal pyramidal patterns is too tight, and
a more general solution must be found.


Conjugated systems
------------------
No detection algorithms are in place to find conjugated systems. Conjugated
systems would require special care in the following places:

- In ranking, the enumeration of Kekule structures is part of a sequence rule
  priority determination algorithm.
- In bond order discretization, aromatic bond orders (around 1.5) may be
  randomly rounded up or down, with no care taken to generate a Kekule structure
  of a contained conjugated system
- In conformer generation, conjugated systems should be flat.
- In molecule comparison, mesomer kekule structures are not recognized as
  equivalent molecules
