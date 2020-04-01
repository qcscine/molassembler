making the distance bounds on 1-3 distances looser does not significantly
improve average error function values for pentagonal bipyramidal or 
square antiprismatic shapes.

including 1-3 distances likely overconstrains these shapes with large
coordination numbers, as can be seen from a degree-of-freedom calculation,
assuming all shapes considered are non-linear:

N = coordination number
number of 1-2 distance constraints = N
number of 1-3 distance constraints = N(N-1)/2
free DOF = 3(N + 1) - 6 - N - N(N-1)/2
         = ½(-N² + 5N - 6)

N         2   3   4   5   6   7   8 
free DOF  0   0  -1  -3  -6 -10 -15

error functions are very likely to increase over such overconstrained systems
whose distance matrices likely do not satisfy triangle inequalities yet.

However, correcting mistakes made in the angle functions does correct the major
increase in average error function value. Although shapes with high
coordination numbers do not reach 0 in optimizations (likely due to triangle
inequality violations), even square pyramidal is now below 0.1 (previously at
4). The generated structures look much better. Getting the refinement to reach
the global minimum is now the target. Introducing metrization may allow this.

Metrization changes nothing. No triangle inequalities are introduced by random
independent choice of distances.
