"""
Insert a carbon in the molecular graph between two atoms that have no formal charge.
"""

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule

from .action_molecular_graph import ActionMolGraph


class InsertCarbonMolGraph(ActionMolGraph):
    """
    Insert a carbon in the molecular graph between two atoms that have no formal charge.
    The new carbon is bonded with two single bonds to the two atoms.
    The initial bond is removed.
    """

    def __init__(
        self,
        molecule: Molecule,
        atom1_idx_to_bond: int,
        bond_to_1: int,
        atom2_idx_to_bond: int,
        bond_to_2: int,
    ) -> None:
        super().__init__(molecule)
        # index of the atom to bond to the new atom
        self.atom1_idx_to_bond: int = atom1_idx_to_bond
        self.bond_to_1: int = bond_to_1
        self.atom2_idx_to_bond: int = atom2_idx_to_bond
        self.bond_to_2: int = bond_to_2

    @override
    def apply_action(self, new_mol_graph: MolecularGraph) -> None:
        # Removing bond between atoms
        new_mol_graph.set_bond(self.atom1_idx_to_bond, self.atom2_idx_to_bond, 0)

        # Adding the atom
        new_mol_graph.add_atom("C")
        new_atom_idx: int = new_mol_graph.nb_atoms - 1

        # Creating a bond from the last inserted atom to the existing ones
        new_mol_graph.set_bond(new_atom_idx, self.atom1_idx_to_bond, self.bond_to_1)
        new_mol_graph.set_bond(new_atom_idx, self.atom2_idx_to_bond, self.bond_to_2)

    def __repr__(self) -> str:
        return (
            f"InsertCarbonMolGraph({self.molecule}, "
            f"{self.atom1_idx_to_bond}, {self.atom2_idx_to_bond})"
        )

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to X to the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        # insertion only possible if the molecule is not at its maximum size
        if mol_graph.nb_atoms >= Molecule.max_heavy_atoms:
            return []

        action_list: list[Action] = []

        charged_or_radical: list[bool] = [
            mol_graph.atom_charged_or_radical(atom)
            for atom in range(mol_graph.nb_atoms)
        ]

        # insert on existing bond of non charged atoms if at least one atom is mutable
        for atom1, atom2 in mol_graph.bonds:

            bond_type = mol_graph.bond_order(atom1, atom2)
            charged_or_radical1 = charged_or_radical[atom1]
            charged_or_radical2 = charged_or_radical[atom2]
            # bond exists between the two atoms
            if bond_type == 0:
                continue
            # if a1 and a2 are charged or radical
            if charged_or_radical1 and charged_or_radical2:
                # C can be inserted only if the C can be connected to the two
                # atoms without changing the bond type
                # limit of 4 bonds for C so bond_type <= 2 for each atom
                if bond_type <= 2:
                    action_list.append(
                        InsertCarbonMolGraph(
                            molecule, atom1, bond_type, atom2, bond_type
                        )
                    )
                continue

            # if only one of the two atoms is charged or radical
            if charged_or_radical1 or charged_or_radical2:
                bond_to_1 = bond_type if charged_or_radical1 else 1
                bond_to_2 = bond_type if charged_or_radical2 else 1
                action_list.append(
                    InsertCarbonMolGraph(molecule, atom1, bond_to_1, atom2, bond_to_2)
                )
                continue
            action_list.append(InsertCarbonMolGraph(molecule, atom1, 1, atom2, 1))
        return action_list
