"""
Insert a carbon in the molecular graph between two atoms that have no formal charge.
"""

import itertools

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, ActionSpace, Molecule


class InsertCarbonMolGraph(Action):
    """
    Insert a carbon in the molecular graph between two atoms that have no formal charge.
    The new carbon is bonded with two single bonds to the two atoms.
    The initial bond is removed.
    """

    def __init__(self, atom1_idx_to_bond: int, atom2_idx_to_bond: int) -> None:
        # index of the atom to bond to the new atom
        self.atom1_idx_to_bond: int = atom1_idx_to_bond
        self.atom2_idx_to_bond: int = atom2_idx_to_bond

    @override
    def apply(self, molecule: Molecule) -> Molecule:
        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph = MolecularGraph(mol_graph.smiles)

        # Removing bond between atoms
        new_mol_graph.set_bond(self.atom1_idx_to_bond, self.atom2_idx_to_bond, 0)

        # Adding the atom
        new_mol_graph.add_atom("C")
        new_atom_idx: int = new_mol_graph.nb_atoms - 1

        # Creating a bond from the last inserted atom to the existing ones
        new_mol_graph.set_bond(new_atom_idx, self.atom1_idx_to_bond, 1)
        new_mol_graph.set_bond(new_atom_idx, self.atom2_idx_to_bond, 1)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Insert Carbon caused an error.")
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return (
            f"InsertCarbonMolGraph({self.atom1_idx_to_bond}, {self.atom2_idx_to_bond})"
        )


class ActionSpaceInsertCarbonMolGraph(ActionSpace):
    """
    List possible actions on molecular graphs to X.
    """

    @override
    def list_actions(self, molecule: Molecule) -> list[Action]:
        """List possible actions to X to the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        # insertion only possible if the molecule is not at its maximum size
        if mol_graph.nb_atoms >= Molecule.max_heavy_atoms:
            return []

        action_list: list[Action] = []

        formal_charge_vector = mol_graph.formal_charge_vector()

        # insert on existing bond of non charged atoms if at least one atom is mutable
        for atom1, atom2 in itertools.combinations(range(mol_graph.nb_atoms), 2):
            if (
                # bond exists between the two atoms
                mol_graph.bond_type_num(atom1, atom2) > 0
                # no formal charge on the two atoms
                and formal_charge_vector[atom1] == 0
                and formal_charge_vector[atom2] == 0
                # at least one atom is mutable
                and (
                    mol_graph.atom_mutability(atom1) or mol_graph.atom_mutability(atom2)
                )
            ):
                action_list.append(InsertCarbonMolGraph(atom1, atom2))
        return action_list
