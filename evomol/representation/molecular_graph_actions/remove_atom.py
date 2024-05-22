"""
Remove an atom from the molecular graph.
"""

import networkx as nx
from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule


class RemoveAtomMolGraph(Action):
    """
    Remove an atom from the molecular graph.
    """

    def __init__(self, molecule: Molecule, atom_idx: int) -> None:
        super().__init__(molecule)
        # index of the atom to remove
        self.atom_idx: int = atom_idx

    @override
    def apply(self) -> Molecule:
        mol_graph: MolecularGraph = self.molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph = MolecularGraph(mol_graph.smiles)

        # Removing the atom
        new_mol_graph.remove_atom(self.atom_idx)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Remove atom caused an error.")
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return f"RemoveAtomMolGraph({self.molecule}, {self.atom_idx})"

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to remove an atom from the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        # impossible to remove the atoms whose removal would create
        # multiple connected components in the molecular graph

        # convert the adjacency matrix to a networkx graph
        adjacency_matrix = nx.from_numpy_array(mol_graph.adjacency_matrix)

        # find articulation points
        # (i.e. removing the atom would create two connected components)
        articulation_points = [False] * mol_graph.nb_atoms
        for atom_idx in nx.articulation_points(adjacency_matrix):
            articulation_points[atom_idx] = True

        # any atom can be removed if it is mutable and
        # if it has only one neighbor or if none of its bonds are bridges
        for i in range(mol_graph.nb_atoms):
            if not articulation_points[i] and mol_graph.atom_mutability(i):
                action_list.append(RemoveAtomMolGraph(molecule, i))

        return action_list
