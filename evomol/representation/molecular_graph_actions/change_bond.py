"""
Change bond for molecular graph representation.
"""

import itertools
from copy import copy

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule


class ChangeBondMolGraph(Action):
    """
    Changing a bond from any type to any type among no bond, single, double, triple.

    """

    avoid_bond_breaking: bool = False

    def __init__(
        self,
        molecule: Molecule,
        atom1: int,
        atom2: int,
        bond_type: int,
    ) -> None:
        super().__init__(molecule)
        self.atom1: int = atom1
        self.atom2: int = atom2
        self.bond_type: int = bond_type

    @override
    def apply(self) -> Molecule:
        mol_graph: MolecularGraph = self.molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph: MolecularGraph = copy(mol_graph)

        new_mol_graph.set_bond(self.atom1, self.atom2, self.bond_type)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Change bond caused an error.")
            raise e

        return Molecule(new_mol_graph.smiles)

    def __repr__(self) -> str:
        return (
            f"ChangeBondMolGraph({self.molecule}, "
            f"{self.atom1}, {self.atom2}, {self.bond_type})"
        )

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to change a bond in the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        implicit_valences = mol_graph.implicit_valences
        bridge_bonds_matrix = mol_graph.bridge_bonds_matrix

        charged_or_radical: list[bool] = [
            mol_graph.atom_charged_or_radical(atom)
            for atom in range(mol_graph.nb_atoms)
        ]

        # for each bond
        for atom1, atom2 in itertools.combinations(range(mol_graph.nb_atoms), 2):
            if charged_or_radical[atom1] or charged_or_radical[atom2]:
                continue

            current_bond: int = mol_graph.bond_type_num(atom1, atom2)

            # the max bond that can be formed between atom1 and atom2
            # is the minimum of the implicit valence of the two atoms
            # plus the current bond
            max_bond: int = (
                min(implicit_valences[atom1], implicit_valences[atom2]) + current_bond
            )

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
                        not bridge_bonds_matrix[atom1][atom2]
                        and not cls.avoid_bond_breaking
                    ):
                        action_list.append(
                            ChangeBondMolGraph(molecule, atom1, atom2, new_bond)
                        )

                # Bond increment
                # Bond can be incremented only if the new bond is less than
                # the maximum bond that can be formed between the two atoms
                elif max_bond >= new_bond:
                    action_list.append(
                        ChangeBondMolGraph(molecule, atom1, atom2, new_bond)
                    )

        return action_list
