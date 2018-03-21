# Strained molecule tests fail

Observations:
1. 9 molecules fail:
   - fenestrane-4
   - bridged-fused-cyclopropanes-ketone
   - strained-db-aromatic-multicycles-1, 2, 3
   - biphenyl-cross-bridged-w-triple-bonds
   - cyclobutadiene
   - fused-cyclopentane-db-at-bridgehead
   - cyclohepta-3,3,7,7-tetramethyl-1-ine


# Sanity tests fail

Observations:
1. Most tests are very close. However, in the cases that I can see, some
   geometries are particularly troublesome in distinguishing. Perhaps that is
   due to not perfect optimized geometries for these examples?
2. Particularly difficult is apparently distinguishing:
   - square planar / seesaw
   - seesaw / seesaw
3. Sanity tests are going to be a heavy computational testing expense if truly
   every single arrangement of the stereocenters are tested (remember
   square-antiprismatic has >5000 for the maximally asymmetric case)
4. In trigonal-prismatic, there is a failure due to exceeding refinement
   failure threshold!

Possible sources of error:

Possible actions:
- Look at output geometries again to see if they are truly optimized into a
  proper minimum
- Make the tests into more of a spot check instead of a full computational
  exercise
- Use real world example molecules, particularly ones that can be checked with
  geometry determination methods (currently only VSEPR) instead of contrived
  unlogical ones
- Chiral centers may need to be stricter than surrounding areas, maybe locally
  different oneTwoVariance parameters may be helpful

Actions taken:

Certainties about source of error:

Conclusions:
- Look at trigonal prismatic again -> The Stereocenter information encapsulation
  through the tetrahedra I set is questionable! Sanity tests fail on that
  symmetry due to high failure rate -> investigate!


# Issue header

Observations:

Possible sources of error:

Possible actions:

Actions taken:

Certainties about source of error:

Fix:

Conclusions:
