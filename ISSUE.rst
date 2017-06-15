Not 100% reliable for more-than-minimal edge length stddevs
-----------------------------------------------------------

Observations:
- I have to discard some roots that are *larger* than the first enclosing 
  circumradius. I shouldn't have to do that, there ought to be seven roots 
  *every* time, and the first of those should be the cyclic pentagon.
- Made a graph of stddev vs. success rate, and past a minimal edge length
  stddev, reliability flounders around 95%.

Possible sources of error:
- Function is incorrect, either due to numerical issues or incorrect calculation

Possible actions:
- Ask paper author to run some of my failing edge sequences through his
  implementation and see if his works. Should yield certainty about where error
  lies, with me - or with the paper.

Actions taken:
- Replace all doubles with long doubles. No effect on failure cases.
- Check whether R version has the same failings. It does. Reinvestigating
  fluctuations close to zero for missing roots, although it makes little sense,
  after all the circumradii there are wildly out of reasonable bounds! Using
  MPFR for smaller circumradii cleans up fluctuations, but reveals no additional
  roots.
- Recheck against paper. Did that, and all coefficients & the general
  calculation scheme seem correct.

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
