"""
SMILES representation of a molecule.
"""

from typing_extensions import override

from evomol.representation.molecule import MoleculeRepresentation


class SMILES(MoleculeRepresentation):
    """Molecule representation using SMILES string."""

    def __init__(self, smiles: str):
        """Initialize the molecule representation using a SMILES string.

        Args:
            smiles (str): SMILES string of the molecule
        """
        super().__init__(smiles)
        self.smiles = smiles

    @classmethod
    def set_generic_action_space(cls) -> None:
        """Set the generic action space for the SMILES representation."""
        cls.set_action_spaces([])

    @override
    def representation(self) -> str:
        return self.smiles

    def __repr__(self) -> str:
        return f"SMILES({self.smiles})"
