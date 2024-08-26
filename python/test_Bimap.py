import scine_molassembler as masm

def test_Bimap_comparison():

    neopentane = masm.io.experimental.from_smiles("CC(C)(C)C")
    methyl = masm.io.experimental.from_smiles("[CH3]")
    matches = masm.subgraphs.complete(methyl, neopentane)
    another_matches = masm.subgraphs.complete(methyl, neopentane)

    assert matches[0] == matches[0]
    assert matches[0] == another_matches[0]
    assert all(i == j for i, j in zip(matches, another_matches))
    assert matches[0] != matches[1]
    assert id(matches[0]) != id(another_matches[0])
