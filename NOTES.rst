TODO
----

- ligand loss algorithm is not complete, neither in dynamic or constexpr
  contexts
- One issue that becomes apparent in the mappings work is that ligand
  equivalencies are not discovered because the symmetry coordinates deviate from
  ideality

  - Replace sq antipr with idealized coordinates (done)
  - Replace equality criteria in algorithm with approximate equality (e.g. 1e-6)
    (ALSO in both ligand gain algorithms!)
