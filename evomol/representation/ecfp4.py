"""
Module for the Extended Connectivity Fingerprint 4 (ECFP4) representation of a molecule.
"""

from evomol.representation.molecule import MoleculeRepresentation


class ECFP4(MoleculeRepresentation):
    """
    Molecule representation using the Extended Connectivity Fingerprint 4 (ECFP4).
    """

    def __init__(self, smiles: str):
        super().__init__(smiles)
        self.smiles = smiles

    def representation(self) -> str:
        raise NotImplementedError
