Numerically imprecise
---------------------

Observations:
- The root-finding idea of scanning the function from 0 outward and using
  TOMS748 on the first crossing fails due to loads of fluctuation at low rhos.
  The function returns values in increments of magnitude 1e-6.

Possible sources of error:
- Summations. Tried adding more kahan summation or replacing it with long
  summation, and although both induced some change, it was not significant in
  any remote fashion.
- Multiplications / Divisions: At most 0.5 ulp, hard to imagine this being
  unstable

Possible actions:
- transition to long doubles everywhere

Actions taken:

Certainties about source of error:

Conclusions:
