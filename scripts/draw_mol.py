import sys
import os


# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from evomol.representation import Molecule, SMILES, MolecularGraph
from evomol.visualization.draw_mol import (
    mol_to_svg,
    draw_multiple_svgs_in_matplotlib,
    draw_in_matplotlib,
)

Molecule.id_representation_class = SMILES
Molecule.representations_class = [MolecularGraph]


if __name__ == "__main__":
    molecules = [Molecule("CCO"), Molecule("CCC"), Molecule("CCN")]

    note_on_atoms = {0: "A", 1: "B", 2: "C"}
    note_on_bonds = {0: "bond A", 1: "bond B"}
    svg_1 = mol_to_svg(
        molecules[0],
        notes_on_atoms=note_on_atoms,
        notes_on_bonds=note_on_bonds,
    )

    draw_in_matplotlib(svg_1)

    svgs = [mol_to_svg(mol) for mol in molecules]
    draw_multiple_svgs_in_matplotlib(svgs)
