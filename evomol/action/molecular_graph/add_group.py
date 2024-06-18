"""
Add a group of atoms to the molecular graph.
"""

from dataclasses import dataclass

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule

from .action_molecular_graph import ActionMolGraph


@dataclass
class Group:
    """
    Group of atoms to add to the molecular graph.
    """

    smiles: str
    nb_atoms: int
    positions_to_link: list[int]


class AddGroupMolGraph(ActionMolGraph):
    """
    Add a group of atoms to the molecular graph.
    """

    groups: list[Group] = [
        Group("C1=CC=CS1", 5, [0]),
        Group("C1=CC=CC=C1", 6, [0]),
        Group("[N+](=O)[O-]", 3, [0]),
        Group("N=[N+]=[N-]", 3, [0]),
        Group("S(=O)(=O)O", 4, [0]),
    ]

    def __init__(
        self,
        molecule: Molecule,
        atom_to_link: int,
        smiles: str,
        added_group_atom: int,
    ) -> None:
        super().__init__(molecule)
        self.atom_to_link: int = atom_to_link
        self.smiles: str = smiles
        self.added_group_atom: int = added_group_atom

    @override
    def apply_action(self, new_mol_graph: MolecularGraph) -> None:
        if new_mol_graph.nb_atoms == 0:
            # If the molecule is empty
            # the new molecular graph is the group
            new_mol_graph.change_smiles(self.smiles)
            return

        # Add the group to the molecular graph
        new_mol_graph.add_group(
            self.atom_to_link,
            self.smiles,
            self.added_group_atom,
        )

    def __repr__(self) -> str:
        return f"AddGroupMolGraph({self.molecule}, {self.smiles}, {self.atom_to_link})"

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to add a group to the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None

        # if the molecule is too big, we cannot add a group
        if mol_graph.nb_atoms >= Molecule.max_heavy_atoms:
            return []

        # if the molecule is empty, we can add any group
        if mol_graph.nb_atoms == 0:
            return [
                AddGroupMolGraph(molecule, 0, group.smiles, position)
                for group in cls.groups
                for position in group.positions_to_link
            ]

        # otherwise, for each atom, if the implicit valence is greater than 0
        # an action can be performed for each possible group
        # and for each possible position to link
        # if the size of the molecule is still below the maximum size
        implicit_valences = mol_graph.implicit_valences

        return [
            AddGroupMolGraph(molecule, atom_idx, group.smiles, position)
            for atom_idx in range(mol_graph.nb_atoms)
            if implicit_valences[atom_idx] > 0
            and not mol_graph.atom_charged_or_radical(atom_idx)
            for group in cls.groups
            if mol_graph.nb_atoms + group.nb_atoms <= Molecule.max_heavy_atoms
            for position in group.positions_to_link
        ]
