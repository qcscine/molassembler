TODO
----

- Maybe it would have been best to originally use DynamicArrays everywhere with
  the maximum symmetry size or reasonable defaults. Consequentially could have
  done everything with much fewer template patterns... If any are exceeded it
  doesn't compile, so you're fine. Plus it could be stored homogeneously, which
  would be a huge plus
- DynamicProperties' generateAllRotations doesn't have the hash optimization
- Merge in geometry_assignment? Would permit constexpr assignment evaluation,
  which would allow constexpr tables of number of assignments given a number of
  adjacent terminal hydrogens (great optimization for RankingTree) for every
  symmetry
- Tests do not cover intra-size or loss mapping situations
- Tests do not point-check correctness of numUnlinkedAssignments, merely
  discrepancies between dynamic and constexpr solutions
