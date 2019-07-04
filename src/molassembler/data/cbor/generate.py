import molassembler

molecules = {
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "biphenyl": "C1=CC=C(C=C1)C2=CC=CC=C2",
    "pyrimidine": "C1=CN=CN=C1",
    "pyridine": "C1=CC=NC=C1",
    "methane": "C",
    "chlorobenzene": "C1=CC=C(C=C1)Cl",
    "phenole": "C1=CC=C(C=C1)O",
    "mesitylene": "CC1=CC(=CC(=C1)C)C",
    "nhc": "C1NC=CN1"
}

for name, smiles in molecules.items():
    mol = molassembler.io.LineNotation.from_canonical_smiles(smiles)
    molassembler.io.write(filename=(name + ".cbor"), molecule=mol)
