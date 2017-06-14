Unreliable
----------

Observations:
- After introducing root-checking, 17/100 fail

Possible sources of error:
- Numerical issues

Possible actions:

Actions taken:
- Replace all doubles with long doubles

Certainties about source of error:

Conclusions:


Numerically imprecise
---------------------

Observations:
- The root-finding idea of scanning the function from 0 outward and using
  TOMS748 on the first crossing fails due to loads of fluctuation at low rhos.
  The function returns values in increments of magnitude 1e-6.
- One consideration that can be made to further bound the sought range for rho
  is that the minimum radius for the cyclic pentagon is the regular pentagon
  constructed using the minimum edge length, and likewise for the maximum. 
  From this, maximum and minimum rhos can be constructed to search for the
  correct root.

Possible sources of error:
- Summations 
- Multiplications / Divisions: At most 0.5 ulp, hard to imagine this being
  unstable

Possible actions:

Actions taken:
- Tried adding more kahan summation or replacing it with long
  summation, and although both induced some change, it was not significant in
  any remote fashion.
- Transition to long doubles everywhere. This works and changes the field
  entirely, yet does not expose additional structure. There is still too much
  fluctuation at low rhos.

Certainties about source of error:

Conclusions:
- Actually not an issue. No crossings at such low rhos are relevant since they
  are ridiculously large cirumradii. Search only between considerations between
  minimum and maximum.
