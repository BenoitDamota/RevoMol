"""
Move a group of atoms in a molecular graph.
"""

from copy import copy

from typing_extensions import override
import rdkit

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule


class MoveGroupMolGraph(Action):
    """
    Move a group of atoms at a different position in the molecular graph.
    """

    def __init__(
        self,
        molecule: Molecule,
        atom_moving: int,
        atom_staying: int,
        atom_to_link: int,
        bond_type: int,
    ) -> None:
        super().__init__(molecule)
        # index of the atoms in the bridge to move
        self.atom_moving: int = atom_moving
        self.atom_staying: int = atom_staying
        # index of the atom to link to the group
        self.atom_to_link: int = atom_to_link
        self.bond_type: int = bond_type

    @override
    def apply(self) -> Molecule:
        mol_graph: MolecularGraph = self.molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph: MolecularGraph = copy(mol_graph)

        # Removing the bond
        new_mol_graph.set_bond(self.atom_moving, self.atom_staying, 0)

        # Setting the new bond
        new_mol_graph.set_bond(self.atom_moving, self.atom_to_link, self.bond_type)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Move group caused an error.")
            print("SMILES before: " + mol_graph.smiles)
            print("SMILES after: " + new_mol_graph.smiles)
            raise e

        return Molecule(new_mol_graph.smiles)

    def __repr__(self) -> str:
        return (
            "MoveGroupMolGraph("
            f"{self.molecule}, {self.atom_moving}, "
            f"{self.atom_staying}, {self.atom_to_link}, {self.bond_type})"
        )

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to move a functional group in the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        nb_atoms = mol_graph.nb_atoms

        charged_or_radical: list[bool] = [
            mol_graph.atom_charged_or_radical(atom)
            for atom in range(mol_graph.nb_atoms)
        ]

        # list bridge bonds where at most one atom is charged or radical
        bridges = [
            (atom1, atom2)
            for atom1, atom2 in mol_graph.bridge_bonds
            if not (charged_or_radical[atom1] and charged_or_radical[atom2])
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

        implicit_valences = mol_graph.implicit_valences

        # for each bridge bond
        for atom1, atom2 in bridges:

            charged_or_radical1 = charged_or_radical[atom1]
            charged_or_radical2 = charged_or_radical[atom2]

            # if one of the atom is charged or radical,
            # the type of bond cannot be changed
            bond_type = (
                1
                if not charged_or_radical1 and not charged_or_radical2
                else mol_graph.bond_type_num(atom1, atom2)
            )

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

            # atom2 as atom_moving and atom1 as atom_staying
            if not charged_or_radical1:
                # try to connect atom2 to the atoms in component1
                # if the implicit valence is enough
                for atom_to_link in component1:
                    if implicit_valences[atom_to_link] >= bond_type:
                        action_list.append(
                            MoveGroupMolGraph(
                                molecule,
                                atom_moving=atom2,
                                atom_staying=atom1,
                                atom_to_link=atom_to_link,
                                bond_type=bond_type,
                            )
                        )

            # atom1 as atom_moving and atom2 as atom_staying
            if not charged_or_radical2:
                # try to connect atom1 to the atoms in component2
                # if the implicit valence is enough
                for atom_to_link in component2:
                    if implicit_valences[atom_to_link] >= bond_type:
                        action_list.append(
                            MoveGroupMolGraph(
                                molecule,
                                atom_moving=atom1,
                                atom_staying=atom2,
                                atom_to_link=atom_to_link,
                                bond_type=bond_type,
                            )
                        )

        return action_list
