"""
Remove a group of connected atoms from the molecular graph.
"""

from copy import copy

import networkx as nx
from typing_extensions import override
import rdkit

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule


def connected_component_after_removing_bond(
    adjacency_matrix: list[list[int]], atom1: int, atom2: int
) -> list[int]:
    """Return the connected component containing atom1 after removing the bond
    between atom1 and atom2.

    Args:
        adjacency_matrix (list[list[int]]): adjacency matrix of the molecular graph
        atom1 (int): first atom index and atom in the connected component return
        atom2 (int): second atom index

    Returns:
        list[int]: list of atom indices in the connected component containing atom1
    """


class RemoveGroupMolGraph(Action):
    """
    Remove a group of connected atoms from the molecular graph.
    """

    remove_only_smallest: bool = False

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
    def apply(self) -> Molecule:
        mol_graph: MolecularGraph = self.molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph: MolecularGraph = copy(mol_graph)

        # Remove the bond between the two bridge atoms
        new_mol_graph.set_bond(self.bridge_atom_to_keep, self.bridge_atom_to_remove, 0)

        # Remove the atoms in the connected component containing
        # the bridge atom to remove in reverse order
        to_remove: list[str] = sorted(
            list(
                nx.node_connected_component(
                    nx.Graph(new_mol_graph.adjacency_matrix), self.bridge_atom_to_remove
                )
            ),
            reverse=True,
        )
        for atom in to_remove:
            new_mol_graph.remove_atom(atom)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Removing group caused an error.")
            print("SMILES before: " + mol_graph.smiles)
            print("SMILES after: " + new_mol_graph.smiles)
            raise e

        return Molecule(new_mol_graph.smiles)

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
            if not charged_or_radical[atom1] and not charged_or_radical[atom2]
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

            any_charged1 = any(charged_or_radical[atom] for atom in component1)
            any_charged2 = any(charged_or_radical[atom] for atom in component2)

            if any_charged1 and any_charged2:
                continue

            if cls.remove_only_smallest:
                if len(component1) <= len(component2) and not any_charged1:
                    action_list.append(
                        RemoveGroupMolGraph(
                            molecule,
                            bridge_atom_to_keep=atom2,
                            bridge_atom_to_remove=atom1,
                        )
                    )
                if len(component2) <= len(component1) and not any_charged2:
                    action_list.append(
                        RemoveGroupMolGraph(
                            molecule,
                            bridge_atom_to_keep=atom1,
                            bridge_atom_to_remove=atom2,
                        )
                    )
            else:
                if not any_charged2:
                    action_list.append(
                        RemoveGroupMolGraph(
                            molecule,
                            bridge_atom_to_keep=atom1,
                            bridge_atom_to_remove=atom2,
                        )
                    )
                if not any_charged1:
                    action_list.append(
                        RemoveGroupMolGraph(
                            molecule,
                            bridge_atom_to_keep=atom2,
                            bridge_atom_to_remove=atom1,
                        )
                    )

        return action_list
