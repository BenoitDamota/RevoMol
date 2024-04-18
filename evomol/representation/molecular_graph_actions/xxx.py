"""
Example of a molecular graph action.
"""

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, ActionSpace, Molecule


class XXXMolGraph(Action):
    """
    XXX to the molecular graph.
    """

    def __init__(self) -> None:
        pass

    @override
    def apply(self, molecule: Molecule) -> Molecule:
        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph = MolecularGraph(mol_graph.smiles)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("XXX caused an error.")
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return f"XXXMolGraph({self})"


class ActionSpaceXXXMolGraph(ActionSpace):
    """
    List possible actions on molecular graphs to XXX.
    """

    def __init__(self) -> None:
        pass

    @override
    def list_actions(self, molecule: Molecule) -> list[Action]:
        """List possible actions to XXX to the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        _ = mol_graph

        return action_list
