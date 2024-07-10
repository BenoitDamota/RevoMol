"""
Remove a group of connected atoms from the molecular graph.
"""

import networkx as nx
import rdkit
from typing_extensions import override

from evomol.action import Action
from evomol.representation import MolecularGraph, Molecule

from .action_molecular_graph import ActionMolGraph


class RemoveGroupMG(ActionMolGraph):
    """
    Remove a group of connected atoms from the molecular graph.
    """

    remove_only_smallest: bool = True
    remove_only_single_bond: bool = False
    remove_charged_or_radical: bool = True

    def __init__(
        self,
        molecule: Molecule,
        bridge_atom_to_keep: int,
        bridge_atom_to_remove: int,
    ) -> None:
        super().__init__(molecule)
        self.bridge_atom_to_keep: int = bridge_atom_to_keep
        self.bridge_atom_to_remove: int = bridge_atom_to_remove

    @override
    def __eq__(self, other: object) -> bool:
        return (
            other.__class__ == RemoveGroupMG
            and self.molecule == other.molecule
            and self.bridge_atom_to_keep == other.bridge_atom_to_keep
            and self.bridge_atom_to_remove == other.bridge_atom_to_remove
        )

    @override
    def __hash__(self) -> int:
        return hash(self.__repr__())

    @override
    def apply_action(self, new_mol_graph: MolecularGraph) -> None:
        # Remove the bond between the two bridge atoms
        new_mol_graph.set_bond(self.bridge_atom_to_keep, self.bridge_atom_to_remove, 0)

        # Remove the atoms in the connected component containing
        # the bridge atom to remove in reverse order
        to_remove: list[int] = sorted(
            list(
                nx.node_connected_component(
                    nx.Graph(new_mol_graph.adjacency_matrix), self.bridge_atom_to_remove
                )
            ),
            reverse=True,
        )
        for atom in to_remove:
            new_mol_graph.remove_atom(atom)

    def __repr__(self) -> str:
        return (
            f"RemoveGroupMolGraph("
            f"{self.molecule}, "
            f"{self.bridge_atom_to_keep}, {self.bridge_atom_to_remove})"
        )

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to remove a group from the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        nb_atoms = mol_graph.nb_atoms
        charged_or_radical: list[bool] = [
            mol_graph.atom_charged_or_radical(atom) for atom in range(nb_atoms)
        ]

        # list bridge bonds where no atom is charged or radical
        bridges = [
            (atom1, atom2)
            for atom1, atom2 in mol_graph.bridge_bonds
            if not charged_or_radical[atom1]
            and not charged_or_radical[atom2]
            and (
                not cls.remove_only_single_bond
                or mol_graph.bond_order(atom1, atom2) == 1
            )
        ]

        # No possible action if no bridge bond
        if not bridges:
            return []

        # Compute the distance matrix
        # to determine the connected components for each bridge bond
        # an atom is in the same component as the atom it is closest to
        # in the bridge bond
        distances: list[list[int]] = rdkit.Chem.rdmolops.GetDistanceMatrix(
            mol_graph.mol
        )

        # for each bond
        for atom1, atom2 in bridges:

            # list connected component for each side
            component1 = []
            component2 = []
            for atom in range(nb_atoms):
                if atom in (atom1, atom2):
                    continue
                if distances[atom1][atom] < distances[atom2][atom]:
                    component1.append(atom)
                else:
                    component2.append(atom)

            can_remove_1: bool = (
                cls.remove_charged_or_radical
                or not any(charged_or_radical[atom] for atom in component1)
            ) and (not cls.remove_only_smallest or len(component1) <= len(component2))

            can_remove_2: bool = (
                cls.remove_charged_or_radical
                or not any(charged_or_radical[atom] for atom in component2)
                and (not cls.remove_only_smallest or len(component2) <= len(component1))
            )

            if can_remove_1:
                action_list.append(
                    RemoveGroupMG(
                        molecule,
                        bridge_atom_to_keep=atom2,
                        bridge_atom_to_remove=atom1,
                    )
                )

            if can_remove_2:
                action_list.append(
                    RemoveGroupMG(
                        molecule,
                        bridge_atom_to_keep=atom1,
                        bridge_atom_to_remove=atom2,
                    )
                )

        return action_list
