"""
Remove an atom from the molecular graph.
"""

from copy import copy

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
        new_mol_graph: MolecularGraph = copy(mol_graph)

        # Removing the atom
        new_mol_graph.remove_atom(self.atom_idx)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Remove atom caused an error.")
            raise e

        return Molecule(new_mol_graph.smiles)

    def __repr__(self) -> str:
        return f"RemoveAtomMolGraph({self.molecule}, {self.atom_idx})"

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to remove an atom from the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        # impossible to remove the atoms whose removal would create
        # multiple connected components in the molecular graph

        # find articulation points
        # (i.e. removing the atom would create two connected components)
        articulation_points = [False] * mol_graph.nb_atoms
        for atom_idx in nx.articulation_points(nx.Graph(mol_graph.adjacency_matrix)):
            articulation_points[atom_idx] = True

        charged_or_radical: list[bool] = [
            mol_graph.atom_charged_or_radical(atom)
            for atom in range(mol_graph.nb_atoms)
        ]

        # any atom can be removed if it is not charged nor radical nor
        # an articulation point and
        # none of its neighbors are radical or charged
        return [
            RemoveAtomMolGraph(molecule, atom)
            for atom in range(mol_graph.nb_atoms)
            if not articulation_points[atom]
            and not charged_or_radical[atom]
            and not any(
                charged_or_radical[neighbor]
                for neighbor in mol_graph.atoms_bonded_to(atom)
            )
        ]
