"""
Remove a group of connected atoms from the molecular graph.
"""

import itertools

import networkx as nx
import numpy as np
from typing_extensions import override

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
    # Remove the bond between atom1 and atom2
    adjacency_matrix[atom1][atom2] = 0
    adjacency_matrix[atom2][atom1] = 0

    return list(
        nx.node_connected_component(
            nx.from_numpy_array(np.array(adjacency_matrix)), atom1
        )
    )


class RemoveGroupMolGraph(Action):
    """
    Remove a group of connected atoms from the molecular graph.
    """

    only_remove_smallest_group: bool = False

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
        new_mol_graph = MolecularGraph(mol_graph.smiles)

        atom_to_remove = connected_component_after_removing_bond(
            mol_graph.adjacency_matrix,
            self.bridge_atom_to_remove,
            self.bridge_atom_to_keep,
        )

        # Remove the atoms in reverse order to avoid changing the index of the
        # atoms to remove
        for atom in sorted(atom_to_remove, reverse=True):
            new_mol_graph.remove_atom(atom)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Removing group caused an error.")
            print("SMILES before: " + str(mol_graph.canonical_smiles))
            print("SMILES after: " + str(new_mol_graph.canonical_smiles))
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

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

        # Extract the bridge matrix
        bridge_bond_matrix = mol_graph.bridge_bonds_matrix()
        formal_charge_vector = mol_graph.formal_charge_vector()

        # for each bond
        for atom1, atom2 in itertools.combinations(range(mol_graph.nb_atoms), 2):
            # The group can be removed only if the bond is a bridge and none
            # of the atoms has a formal charge and at least one atom is mutable
            if not (
                bridge_bond_matrix[atom1][atom2]
                and formal_charge_vector[atom1] == 0
                and formal_charge_vector[atom2] == 0
                and (
                    mol_graph.atom_mutability(atom1) or mol_graph.atom_mutability(atom2)
                )
            ):
                continue

            if cls.only_remove_smallest_group:
                # extract the indices of both connected components if the current
                # bond was removed
                connected_component_1 = connected_component_after_removing_bond(
                    mol_graph.adjacency_matrix, atom1, atom2
                )
                connected_component_2 = connected_component_after_removing_bond(
                    mol_graph.adjacency_matrix, atom2, atom1
                )
                if len(connected_component_1) <= len(connected_component_2):
                    action_list.append(
                        RemoveGroupMolGraph(
                            molecule,
                            bridge_atom_to_keep=atom2,
                            bridge_atom_to_remove=atom1,
                        )
                    )
                if len(connected_component_2) <= len(connected_component_1):
                    action_list.append(
                        RemoveGroupMolGraph(
                            molecule,
                            bridge_atom_to_keep=atom1,
                            bridge_atom_to_remove=atom2,
                        )
                    )
            else:
                action_list.append(
                    RemoveGroupMolGraph(
                        molecule,
                        bridge_atom_to_keep=atom1,
                        bridge_atom_to_remove=atom2,
                    )
                )
                action_list.append(
                    RemoveGroupMolGraph(
                        molecule,
                        bridge_atom_to_keep=atom2,
                        bridge_atom_to_remove=atom1,
                    )
                )

        return action_list
