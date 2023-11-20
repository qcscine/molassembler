import scine_molassembler as masm

def test_Cycles():
    if not masm.io.LineNotation.enabled:
        return
    
    # GitHub Issue 13 had cycle edges in triplicate due to a bad mix of Python
    # list semantics and C++ iterators returning references
    smiles = "[O]N(O1)O[Bi]123(ON([O])O2)ON([O])O3"
    mol = masm.io.experimental.from_smiles(smiles)
    cycles = [c for c in mol.graph.cycles]
    edges = sorted([edge for cycle in cycles for edge in cycle])
    # There should be no duplicate edges for this particular molecule
    assert len(edges) == len(set(edges))
