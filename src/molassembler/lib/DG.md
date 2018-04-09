# DG procedure details

- If the specified molecule has stereocenters with exactly zero permissible
  stereopermutations, abort immediately. Since this may occur naturally in
  intermediate stages of graph modification, this is generally permitted.
  Generating a 3D structure for a Molecule that has such a stereocenter is
  however impossible.

# Failure case handling improvement
- If a distance bound violates a triangle inequality, you can recover which
  bound was responsible and could hypothetically loosen any bounds on those
  atoms and re-try before failing out
  This would mean that generateDistanceBounds and generateDistances must return
  a variant, either the success type or a type that contains information on
  where bounds must be loosened
- If propagate fails, then you can try a general loosening several times before
  quitting
