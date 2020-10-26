Directed conformer generator
============================

Making sure you have explored all dihedral angles of a molecule is tricky.
Ideally, a library would approximate the local dihedral energy distribution to
correctly guess the minima and generate those conformers for you. Unfortunately,
Molassembler is short on information needed for such approximations to be
correct. For one, Molassembler avoids requiring full correctness of graph bond
orders to keep errors in floating-point bond order discretization from
propagating. Additionally, Molassembler doesn't ask you to specify the
overall charge of your molecule. The phenomenological approach to local
molecular shape doesn't require it. It is impossible to assign formal charges
purely based on local molecular shapes without knowing overall charge or the
guarantee that the formal bond orders of the graph are correct.

So what can Molassembler offer instead? Given that Molassembler quite explicitly
states that its generated conformers are guesses to local energy minima on the
potential energy surface, what Molassembler can try to ensure is that
minimizations of generated conformers minimize into all local minima that exist.

To that end, we suggest you carefully read the details of
:class:`~scine_molassembler.BondStereopermutator.Alignment` and consider
deduplicating energy minimized conformer guesses with the
:class:`~scine_molassembler.DirectedConformerGenerator.Relabeler`.

.. autoclass:: scine_molassembler.DirectedConformerGenerator
