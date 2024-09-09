"""
This is an example of a new representation class that could be added to the package.
"""

import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from typing_extensions import override

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# pylint: disable=wrong-import-position, import-error

from evomol.representation import SMILES, Molecule, MoleculeRepresentation


class ECFP4(MoleculeRepresentation):
    """
    Molecule representation using the Extended Connectivity Fingerprint 4 (ECFP4).
    """

    def __init__(self, smiles: str):
        super().__init__(smiles)

        mol = Chem.MolFromSmiles(smiles)

        # create a Morgan fingerprints generator
        self.fingerprints_generator = AllChem.GetMorganGenerator(radius=2)
        fingerprints = self.fingerprints_generator.GetSparseCountFingerprint(
            mol
        ).GetNonzeroElements()

        self.fingerprints = set(fingerprints.keys())

    @override
    def representation(self) -> str:
        return str(self.fingerprints)

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, ECFP4):
            return False
        return self.fingerprints == value.fingerprints


def main() -> None:
    """
    Main function to test the ECFP4 representation.
    """

    # initialize the molecule representation
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [ECFP4]
    Molecule.max_heavy_atoms = 6

    Molecule.accepted_atoms = ["C", "O", "N", "S"]

    smiles = "CCO"

    molecule = Molecule(smiles)

    print(molecule.get_representation(ECFP4).representation())


if __name__ == "__main__":
    main()
