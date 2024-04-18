"""
Change bond for molecular graph representation.
"""

import itertools

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, ActionSpace, Molecule


class ChangeBondMolGraph(Action):
    """
    Changing a bond from any type to any type among no bond, single, double, triple.

    """

    def __init__(self, atom1: int, atom2: int, bond_type: int) -> None:
        self.atom1: int = atom1
        self.atom2: int = atom2
        self.bond_type: int = bond_type

    @override
    def apply(self, molecule: Molecule) -> Molecule:
        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph = MolecularGraph(mol_graph.smiles)

        new_mol_graph.set_bond(self.atom1, self.atom2, self.bond_type)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Change bond caused an error.")
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return f"ChangeBondMolGraph({self.atom1}, {self.atom2}, {self.bond_type})"


class ActionSpaceChangeBondMolGraph(ActionSpace):
    """
    List possible actions on molecular graphs to change a bond between two atoms.
    """

    def __init__(self, prevent_removing_creating_bonds: bool = False):
        # whether to prevent the change of bond from type>= 1 to type 0
        # (=breaking the bond)
        self.prevent_removing_creating_bonds: bool = prevent_removing_creating_bonds

    @override
    def list_actions(self, molecule: Molecule) -> list[Action]:
        """List possible actions to change a bond in the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        free_electons_vect = mol_graph.free_electrons_vector()
        formal_charge_vect = mol_graph.formal_charge_vector()
        bridge_matrix = mol_graph.bridge_bonds_matrix()

        # for each bond
        for atom1, atom2 in itertools.combinations(range(mol_graph.nb_atoms), 2):
            current_bond = mol_graph.bond_type_num(atom1, atom2)

            # for each bond type
            for bond_to_form in range(4):
                delta_bond = bond_to_form - current_bond

                formal_charge_ok = (
                    formal_charge_vect[atom1] == 0 and formal_charge_vect[atom2] == 0
                )
                mutability_ok = mol_graph.atom_mutability(
                    atom1
                ) or mol_graph.atom_mutability(atom2)

                if not mutability_ok or not formal_charge_ok:
                    continue

                # Bond decrement
                # only bond that are not bridges can be completely removed.
                # Bond involving atoms with formal charges cannot be changed.
                # Bonds can be changed only if at least one of the atoms is
                # mutable
                if delta_bond < 0:
                    # Checking if the prevent breaking bonds constraint is
                    # respected if set
                    prevent_breaking_bonds_constraint_respected = (
                        not self.prevent_removing_creating_bonds or bond_to_form > 0
                    )

                    if (
                        not bridge_matrix[atom1][atom2] or bond_to_form > 0
                    ) and prevent_breaking_bonds_constraint_respected:
                        action_list.append(
                            ChangeBondMolGraph(atom1, atom2, bond_to_form)
                        )

                # Bond increment
                # Bond can be incremented of delta if each atom involved has at
                # least delta free electrons
                # Bonds involving atoms with formal charges cannot be changed
                elif delta_bond > 0:
                    # Checking if the prevent breaking bonds constraint is
                    # respected if set
                    prevent_breaking_bonds_constraint_respected = (
                        not self.prevent_removing_creating_bonds or current_bond > 0
                    )

                    if (
                        min(free_electons_vect[atom1], free_electons_vect[atom2])
                        >= delta_bond
                        and prevent_breaking_bonds_constraint_respected
                    ):
                        action_list.append(
                            ChangeBondMolGraph(atom1, atom2, bond_to_form)
                        )

        return action_list
