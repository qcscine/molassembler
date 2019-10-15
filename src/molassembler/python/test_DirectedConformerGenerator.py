import pytest

import molassembler as masm
import numpy


def test_DirectedConformerGenerator():
    if not masm.io.LineNotation.enabled:
        return

    arginine_smiles = "C(CC(C(=O)O)N)CN=C(N)N"
    arginine = masm.io.LineNotation.from_isomeric_smiles(arginine_smiles)

    # Create the generator, a decision list, and a conformer
    generator = masm.DirectedConformerGenerator(arginine)
    decision = generator.generate_decision_list()
    conf = None
    while not isinstance(conf, numpy.ndarray):
        conf = generator.generate_conformation(decision)

    # Reinterpret the conformer and make sure it matches
    reinterpreted_decisions = generator.get_decision_list(conf)
    assert reinterpreted_decisions == decision

    # Repeat once
    another_decision = generator.generate_decision_list()
    assert another_decision != decision
    conf = None
    while not isinstance(conf, numpy.ndarray):
        conf = generator.generate_conformation(another_decision)

    # Reinterpret and make sure it matches
    another_reinterpreted_decision = generator.get_decision_list(conf)
    assert another_reinterpreted_decision == another_decision
    assert another_reinterpreted_decision != decision
