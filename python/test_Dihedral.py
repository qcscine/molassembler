import os
from typing import Dict, List, Union
import pytest
import scine_molassembler as masm
import scine_utilities as utils


def get_dihedral_atom_list(filename: str) -> List[int]:
    atoms, bo_collection = utils.io.read(filename)

    masm_results = masm.interpret.molecules(atoms, bo_collection, set(), {}, masm.interpret.BondDiscretization.Binary)

    # Split the atom collection into separate collections for each molecule
    positions = masm_results.component_map.apply(atoms)

    for i, m in enumerate(masm_results.molecules):
        ordering = m.canonicalize(masm.AtomEnvironmentComponents.All)
        # Reorder the molecule-specific atom collection according to the
        # canonicalization reordering
        positions[i] = m.apply_canonicalization_map(ordering, positions[i])

    for i, (m, p) in enumerate(zip(masm_results.molecules, positions)):
        alignment = masm.BondStereopermutator.Alignment.EclipsedAndStaggered
        generator = masm.DirectedConformerGenerator(m, alignment)
        relabeler = generator.relabeler()
        dihedrals = relabeler.add(p.positions)
        dihedral_atoms = []
        for j, d in enumerate(dihedrals):
            dihedral_atom_tuple = (
                relabeler.sequences[j].i_set,
                relabeler.sequences[j].j,
                relabeler.sequences[j].k,
                relabeler.sequences[j].l_set)
            dihedral_atoms.append(dihedral_atom_tuple)

    return dihedral_atoms


def test_neighboring_atoms_are_identical():
    assert get_dihedral_atom_list(
        os.path.join(
            os.path.dirname(
                os.path.realpath(__file__)),
            '../test/data/dihedral/1.mol')) == get_dihedral_atom_list(
        os.path.join(
            os.path.dirname(
                os.path.realpath(__file__)),
            '../test/data/dihedral/2.mol'))
