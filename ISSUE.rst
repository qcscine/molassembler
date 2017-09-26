Compile-time algorithm cost
---------------------------

Observations:
- Pre-computing symmetry mappings for adjacent pairs of symmetries becomes
  excessively expensive for symmetries of size N > 5
- Separating compilation of compile-time algorithm and execution of compile-time
  algorithm has been difficult and left me feeling inconclusive about how it
  works. I believe constexpr functions are compiled first, then optimized (if
  specified), then run.
- Creating pointers to function templates does not incur the heavy compilation
  cost, but evaluating them does. This makes sense, since as long as a function
  is not called, it does not need to exist in any form besides as an address.
- When compiling a program that does not explicitly force the compiler to *run*
  the constexpr functions, merely to *compile* them, a heavy compilation cost is
  incurred (at least 15 min). Maybe in this odd case of
  pointer-to-function-template handling the algorithms are evaluated to check
  correctness in all cases at compile time and invoking the same cost?
- Compiling the tests (at O2), which necessarily also includes instantiating and
  hence compiling the constexpr algorithm, takes 40 seconds
- Merely running the tests (at O2), which involves calling the dynamic algorithm
  and the constexpr algorithm at runtime, takes roughly 3 minutes.
- The constexpr results are handled very similarly between the tests and the
  pre-computation, with the tests being more complicated in principle.
- When compilation is separated into instantiation and evaluation of all mapping
  algorithms by first instantiating all into pointer-to-functions and then
  separately evaluating, templight++ indicates that the most time is spent
  compiling the algorithm rather than evaluating it. This may be misleading,
  though, since the presence of the evaluating code may trigger a more complete
  compilation of the algorithm itself.

Key questions:
- If the tests can compile the constexpr algorithm quickly without evaluation,
  why should the compilation be the cause of the expended time when the program
  forces pre-computation?
- If the evaluation of the tests, which includes the constexpr algorithm
  evaluation *at runtime* is quick, why should the evaluation *at compile-time*
  be slow?

Consequences:
- Assuming the compilation and optimization of the constexpr algorithm is
  equally fast in both use cases (which seems reasonable), the vast difference
  in execution of the constexpr algorithm must be caused by the different
  contexts: Fast at run-time, slow at compile-time.
- It may be caused by a compiler bug, in which case finding out where the issue
  lies is largely impossible

New observations:
- It must be a mixture of both effects. Execution at compile-time is slower, but
  algorithmic improvements at low levels has marked impact on compile time too!
  Focussing on often-repeated code learned from run-time profiling allows me to
  compile in the results for symmetries of size 6 in 15 mins.
