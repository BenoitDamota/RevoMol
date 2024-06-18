"""
Example of a molecular graph action.
"""

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule

from .action_molecular_graph import ActionMolGraph


class XXXMolGraph(ActionMolGraph):
    """
    XXX to the molecular graph.
    """

    # def __init__(self, molecule: Molecule) -> None:
    #     super().__init__(molecule)

    @override
    def apply_action(self, new_mol_graph: MolecularGraph) -> None:
        _ = new_mol_graph

    def __repr__(self) -> str:
        return f"XXXMolGraph({self.molecule})"

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to XXX to the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None

        action_list: list[Action] = []

        _ = mol_graph

        return action_list
