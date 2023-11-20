__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import scine_molassembler as masm

molecules = {
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "biphenyl": "C1=CC=C(C=C1)C2=CC=CC=C2",
    "pyrimidine": "C1=CN=CN=C1",
    "pyridine": "C1=CC=NC=C1",
    "methane": "C",
    "chlorobenzene": "C1=CC=C(C=C1)Cl",
    "phenole": "C1=CC=C(C=C1)O",
    "mesitylene": "CC1=CC(=CC(=C1)C)C",
    "nhc": "C1NC=CN1",
    "multidentate_ligand": "CC(C)(C)P(CC1=NC(=CC=C1)CP(C(C)(C)C)C(C)(C)C)C(C)(C)C",
    "haptic_ligand": "CC(C)(C)N[Si](C)(C)C1=C[C+]C=C1"
}

for name, smiles in molecules.items():
    mol = masm.io.LineNotation.from_canonical_smiles(smiles)
    masm.io.write(filename=(name + ".cbor"), molecule=mol)
