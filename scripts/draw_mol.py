import os
import sys

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from evomol.representation import SMILES, MolecularGraph, Molecule
from evomol.visualization.draw_mol import (
    draw_in_matplotlib,
    draw_multiple_svgs_in_matplotlib,
    mol_to_svg,
)

Molecule.id_representation_class = SMILES
Molecule.representations_class = [MolecularGraph]


if __name__ == "__main__":
    molecules = [
        Molecule("CCO"),
        Molecule("CCC"),
        Molecule("CCN"),
        Molecule("CCCCCC"),
        Molecule("C1=CSC(=C2SC=CS2)S1"),
        Molecule("N1=S=NC2=C1N=S=N2"),
        Molecule(
            "C1=CC(=CC=C1C2=C3C=CC(=C(C4=NC(=C(C5=CC=C(N5)C(=C6C=CC2=N6)C7=CC="
            "C(C=C7)C(=O)O)C8=CC=C(C=C8)C(=O)O)C=C4)C9=CC=C(C=C9)C(=O)O)N3)C(=O)O"
        ),
    ]

    note_on_atoms = {0: "A", 1: "B", 2: "C"}
    note_on_bonds = {0: "bond A", 1: "bond B"}
    svg_1 = mol_to_svg(
        molecules[0],
        notes_on_atoms=note_on_atoms,
        notes_on_bonds=note_on_bonds,
    )

    draw_in_matplotlib(
        svg_1,
        show=False,
        save_to_path="output/visualization/CCO.png",
    )

    svgs = [mol_to_svg(mol) for mol in molecules]
    draw_multiple_svgs_in_matplotlib(
        svgs,
        show=False,
        save_to_path="output/visualization/plot.png",
    )
