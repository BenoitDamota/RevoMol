"""
Count the number of unknown GCF and ECFP for the molecule given in input.

example:
    python scripts/molecule_validity.py "O=S(=O)(O)c1ccccc1c1ccsn1"
"""

import os
import sys

import typer

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol.representation import MolecularGraph, Molecule


def check_molecule(smiles: str) -> None:
    """Check the validity of a molecule.

    Args:
        smiles (str): SMILES representation of the molecule.
    """
    dp.setup_default_parameters()

    evaluations = dp.chembl_and_chembl_zinc_evaluations()

    mol = Molecule(smiles)

    print(
        "canonical no aromatic SMILES:",
        mol.get_representation(MolecularGraph).canonical_smiles,
    )
    print(
        "canonical aromatic SMILES   :",
        mol.get_representation(MolecularGraph).aromatic_canonical_smiles,
    )

    for eval_ in evaluations:
        print(f"{eval_.name:33} : {eval_.evaluate(mol)}")


if __name__ == "__main__":

    typer.run(check_molecule)

    # smiles = "O=S(=O)(O)c1ccccc1c1ccsn1"

    # check_molecule(smiles)
