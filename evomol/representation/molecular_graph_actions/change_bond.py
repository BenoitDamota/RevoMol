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

    def __init__(self, avoid_break_bond: bool = False):
        # whether to prevent the change of bond from type>= 1 to type 0
        # (=breaking the bond)
        self.avoid_break_bond: bool = avoid_break_bond

    @override
    def list_actions(self, molecule: Molecule) -> list[Action]:
        """List possible actions to change a bond in the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        free_electrons = mol_graph.free_electrons_vector()
        formal_charges = mol_graph.formal_charge_vector()
        bridge_bond_matrix = mol_graph.bridge_bonds_matrix()

        # for each bond
        for atom1, atom2 in itertools.combinations(range(mol_graph.nb_atoms), 2):
            current_bond = mol_graph.bond_type_num(atom1, atom2)

            mutability_ok = mol_graph.atom_mutability(
                atom1
            ) or mol_graph.atom_mutability(atom2)

            # Bond involving atoms with formal charges cannot be changed
            formal_charge_ok = formal_charges[atom1] == 0 and formal_charges[atom2] == 0

            if not mutability_ok or not formal_charge_ok:
                continue

            # for each bond type
            for new_bond in (0, 1, 2, 3):
                if current_bond == new_bond:
                    continue

                # Bond decrement
                # only bond that are not bridges can be completely removed.
                # Bonds can be changed only if at least one of the atoms is
                # mutable
                if new_bond < current_bond:
                    if new_bond > 0 or (
                        not bridge_bond_matrix[atom1][atom2]
                        and not self.avoid_break_bond
                    ):
                        action_list.append(ChangeBondMolGraph(atom1, atom2, new_bond))

                # Bond increment
                # delta = new_bond - current_bond
                # Bond can be incremented of delta if each atom involved has at
                # least delta free electrons
                else:
                    delta = new_bond - current_bond
                    if min(free_electrons[atom1], free_electrons[atom2]) >= delta:
                        action_list.append(ChangeBondMolGraph(atom1, atom2, new_bond))

        return action_list
