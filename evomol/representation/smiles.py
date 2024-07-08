"""
SMILES representation of a molecule.
"""

from typing_extensions import override

from evomol.representation.molecule import MoleculeRepresentation


def link(uri: str, label: str = "") -> str:
    """Create a hyperlink in the terminal."""
    if not label:
        label = uri
    parameters = ""

    # OSC 8 ; params ; URI ST <name> OSC 8 ;; ST
    escape_mask = "\033]8;{};{}\033\\{}\033]8;;\033\\"

    return escape_mask.format(parameters, uri, label)


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
        cls.set_action_space([])

    @override
    def representation(self) -> str:
        return self.smiles if self.smiles else "Empty SMILES"
        # return (
        #     f"{
        # link('http://hulab.rxnfinder.org/smi2img/' + self.smiles,
        # self.smiles)}"
        #     if self.smiles
        #     else "Empty SMILES"
        # )

    def __repr__(self) -> str:
        return '"' + self.smiles + '"'
        # return f'SMILES("{self.smiles}")'
        # return f'SMILES("{self.representation()}")'

    def __eq__(self, value: object) -> bool:
        return isinstance(value, SMILES) and self.smiles == value.smiles
