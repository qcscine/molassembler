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
