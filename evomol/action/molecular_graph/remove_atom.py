"""
The action RemoveAtom consist in removing an atom from the molecular graph.
The atom to remove is selected among the atoms that don't represent a bridge
bond (i.e. their removal would not create multiple connected components in the
molecular graph), that are not charged nor radical, and that none of their
neighbors are radical or charged.
"""

import networkx as nx
from typing_extensions import override

from evomol.action import Action
from evomol.representation import MolecularGraph, Molecule

from .action_molecular_graph import ActionMolGraph


class RemoveAtomMG(ActionMolGraph):
    """
    Remove an atom from the molecular graph.
    """

    def __init__(self, molecule: Molecule, atom_idx: int) -> None:
        """Initialize the action to remove an atom from the molecular graph.

        Args:
            molecule (Molecule): Molecule to modify.
            atom_idx (int): Index of the atom to remove.
        """
        super().__init__(molecule)
        # index of the atom to remove
        self.atom_idx: int = atom_idx

    @override
    def apply_action(self, new_mol_graph: MolecularGraph) -> None:
        """Remove the atom from the molecular graph.

        Args:
            new_mol_graph (MolecularGraph): Molecular graph to modify.
        """
        new_mol_graph.remove_atom(self.atom_idx)

    def __repr__(self) -> str:
        return f"RemoveAtomMolGraph({self.molecule}, {self.atom_idx})"

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to remove an atom from the molecular graph.

        Args:
            molecule (Molecule): Molecule to modify.

        Returns:
            list[Action]: List of possible actions to remove an atom.
        """

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
            RemoveAtomMG(molecule, atom)
            for atom in range(mol_graph.nb_atoms)
            if not articulation_points[atom]
            and not charged_or_radical[atom]
            and not any(
                charged_or_radical[neighbor]
                for neighbor in mol_graph.atoms_bonded_to(atom)
            )
        ]
