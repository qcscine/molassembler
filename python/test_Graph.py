import scine_molassembler as masm

def test_graph_contains():

    m = masm.Molecule()
    assert masm.BondIndex(0, 1) in m.graph
    assert masm.BondIndex(1, 2) not in m.graph
