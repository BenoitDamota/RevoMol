"""
Move a functional group in a molecular graph.
"""

import itertools

import networkx as nx
from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, ActionSpace, Molecule


class MoveFunctionalGroupMolGraph(Action):
    """
    Move a functional group at a different position in the molecular graph.
    """

    def __init__(
        self,
        bridge_atom1: int,
        bridge_atom2: int,
        new_link_idx: int,
    ) -> None:
        # index of the atoms in the bridge to move
        self.bridge_atom1: int = bridge_atom1
        self.bridge_atom2: int = bridge_atom2
        # index of the atom to link to the group
        self.new_link_idx: int = new_link_idx

    @override
    def apply(self, molecule: Molecule) -> Molecule:
        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph = MolecularGraph(mol_graph.smiles)

        # compute imputed adjacency matrix
        imputed_adjacency_matrix = mol_graph.adjacency_matrix
        imputed_adjacency_matrix[self.bridge_atom1][self.bridge_atom2] = False
        imputed_adjacency_matrix[self.bridge_atom2][self.bridge_atom1] = False

        # extract the connected components with the bond removed
        connected_components = list(
            nx.connected_components(nx.from_numpy_array(imputed_adjacency_matrix))
        )

        # Select the atom from the initial bond to be bonded to the new atom
        atom_to_bond: int = self.bridge_atom1

        bridge1_in_first_component = self.bridge_atom1 in connected_components[0]
        new_link_in_first_component = self.new_link_idx in connected_components[0]

        if bridge1_in_first_component == new_link_in_first_component:
            atom_to_bond = self.bridge_atom2

        # Extracting the bond type to be removed
        bond_type: int = new_mol_graph.bond_type_num(
            self.bridge_atom1, self.bridge_atom2
        )

        # Removing the bond
        new_mol_graph.set_bond(self.bridge_atom1, self.bridge_atom2, 0)

        # Setting the new bond
        new_mol_graph.set_bond(atom_to_bond, self.new_link_idx, bond_type)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Move group caused an error.")
            print("SMILES before: " + mol_graph.canonical_smiles)
            print("SMILES after: " + new_mol_graph.canonical_smiles)
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return (
            "MoveFunctionalGroupMolGraph("
            f"{self.bridge_atom1}, {self.bridge_atom2}, {self.new_link_idx})"
        )


class ActionSpaceMoveFunctionalGroupMolGraph(ActionSpace):
    """
    List possible actions on molecular graphs to move a functional group.
    """

    @override
    def list_actions(self, molecule: Molecule) -> list[Action]:
        """List possible actions to move a functional group in the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        # Extraction of the bridge matrix and free electrons vector
        bridge_bond_matrix = mol_graph.bridge_bonds_matrix()
        free_electrons_vector = mol_graph.free_electrons_vector()
        formal_charge_vector = mol_graph.formal_charge_vector()

        # for each bond
        for atom1, atom2 in itertools.combinations(range(mol_graph.nb_atoms), 2):
            # functional group can be moved only if the bond is a bridge and
            # none of the atoms has a formal charge and at least one atom is mutable
            formal_charge_ok = (
                formal_charge_vector[atom1] == 0 and formal_charge_vector[atom2] == 0
            )
            mutability_ok = mol_graph.atom_mutability(
                atom1
            ) or mol_graph.atom_mutability(atom2)
            if not (
                bridge_bond_matrix[atom1][atom2] and formal_charge_ok and mutability_ok
            ):
                continue
            # Extracting the current bond type
            bond_type_num = mol_graph.bond_type_num(atom1, atom2)

            # for all other atoms
            for k in range(mol_graph.nb_atoms):
                if k in {atom1, atom2}:
                    continue
                # the functional group can be moved if the current atom has
                # enough electrons left and has no formal charge
                if (
                    free_electrons_vector[k] >= bond_type_num
                    and formal_charge_vector[k] == 0
                ):
                    action_list.append(MoveFunctionalGroupMolGraph(atom1, atom2, k))

        return action_list
