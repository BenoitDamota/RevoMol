"""
Cut an atom from the molecular graph.
"""

from typing_extensions import override

from evomol.action import Action
from evomol.representation import MolecularGraph, Molecule

from .action_molecular_graph import ActionMolGraph


class CutAtomMolGraph(ActionMolGraph):
    """
    Cut an atom from the molecular graph.
    """

    def __init__(
        self,
        molecule: Molecule,
        atom_to_cut: int,
        atom1_to_bond: int,
        atom2_to_bond: int,
        bond_type: int,
    ):
        super().__init__(molecule)
        # index of the atom to cut
        self.atom_to_cut: int = atom_to_cut
        # index of the atom to bond to the new atom
        self.atom1_to_bond: int = atom1_to_bond
        self.atom2_to_bond: int = atom2_to_bond
        self.bond_type: int = bond_type

    @override
    def apply_action(self, new_mol_graph: MolecularGraph) -> None:
        # remove bonds
        new_mol_graph.set_bond(self.atom_to_cut, self.atom1_to_bond, 0)
        new_mol_graph.set_bond(self.atom_to_cut, self.atom2_to_bond, 0)

        # create new bond
        new_mol_graph.set_bond(self.atom1_to_bond, self.atom2_to_bond, self.bond_type)

        new_mol_graph.remove_atom(self.atom_to_cut)

    def __repr__(self) -> str:
        return (
            "CutAtomMolGraph("
            f"{self.molecule}, {self.atom_to_cut}, "
            f"{self.atom1_to_bond}, {self.atom2_to_bond}, {self.bond_type})"
        )

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to cut an atom from the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        implicit_valences = mol_graph.implicit_valences
        charged_or_radical: list[bool] = [
            mol_graph.atom_charged_or_radical(atom)
            for atom in range(mol_graph.nb_atoms)
        ]

        for atom_to_cut in range(mol_graph.nb_atoms):
            # we cannot cut a charged or radical atom
            if charged_or_radical[atom_to_cut]:
                continue

            # if there is exactly two atoms bonded to the atom to cut
            bonds_to = mol_graph.atoms_bonded_to(atom_to_cut)
            if len(bonds_to) != 2:
                continue

            atom_1 = bonds_to[0]
            atom_2 = bonds_to[1]

            # atom_1 to atom_to_cut
            bond_1 = mol_graph.bond_order(atom_to_cut, atom_1)

            # atom_2 to atom_to_cut
            bond_2 = mol_graph.bond_order(atom_to_cut, atom_2)

            # atom_1 to atom_2
            bond_12 = mol_graph.bond_order(atom_1, atom_2)

            charged_or_radical1 = charged_or_radical[atom_1]
            charged_or_radical2 = charged_or_radical[atom_2]

            if charged_or_radical1 and charged_or_radical2:
                if bond_1 != bond_2:
                    continue
                action_list.append(
                    CutAtomMolGraph(
                        molecule,
                        atom_to_cut,
                        atom_1,
                        atom_2,
                        bond_1 + bond_12,
                    )
                )
                continue

            max_valence_1 = implicit_valences[atom_1] + bond_1
            max_valence_2 = implicit_valences[atom_2] + bond_2

            if charged_or_radical1 and max_valence_2 >= bond_1:
                action_list.append(
                    CutAtomMolGraph(
                        molecule,
                        atom_to_cut,
                        atom_1,
                        atom_2,
                        bond_1 + bond_12,
                    )
                )
                continue

            if charged_or_radical2 and max_valence_1 >= bond_2:
                action_list.append(
                    CutAtomMolGraph(
                        molecule,
                        atom_to_cut,
                        atom_1,
                        atom_2,
                        bond_2 + bond_12,
                    )
                )
                continue

            if charged_or_radical1 or charged_or_radical2:
                continue

            action_list.append(
                CutAtomMolGraph(
                    molecule,
                    atom_to_cut,
                    atom_1,
                    atom_2,
                    bond_12 + min(bond_1, bond_2),
                )
            )

        return action_list
