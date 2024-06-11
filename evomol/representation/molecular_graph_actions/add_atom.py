"""
Add an atom to a molecular graph.
"""

from copy import copy

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule


class AddAtomMolGraph(Action):
    """
    Add an atom to the molecular graph.
    """

    def __init__(self, molecule: Molecule, index_atom: int, atom_type: str) -> None:
        super().__init__(molecule)
        # index of the atom to bond to the new atom
        self.index_atom: int = index_atom
        # type of the new atom
        self.atom_type: str = atom_type

    @override
    def apply(self) -> Molecule:
        mol_graph: MolecularGraph = self.molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph: MolecularGraph = copy(mol_graph)

        # Adding the atom
        new_mol_graph.add_atom(self.atom_type)

        if new_mol_graph.nb_atoms > 1:
            # Creating a bond from the last inserted atom to the existing one
            new_mol_graph.set_bond(new_mol_graph.nb_atoms - 1, self.index_atom, 1)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Addition caused an error.")
            raise e

        return Molecule(new_mol_graph.smiles)

    def __repr__(self) -> str:
        return f"AddAtomMolGraph({self.molecule}, {self.index_atom}, {self.atom_type})"

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to add an atom to the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None

        # if the molecule is too big, we cannot add an atom
        if mol_graph.nb_atoms >= Molecule.max_heavy_atoms:
            return []

        # if the molecule is empty, we can add any atom
        if mol_graph.nb_atoms == 0:
            return [
                AddAtomMolGraph(molecule, 0, atom_type)
                for atom_type in Molecule.accepted_atoms
            ]

        # otherwise, for each atom, if the implicit valence is greater than 0
        # an action can be performed for each possible atom type

        implicit_valences = mol_graph.implicit_valences

        return [
            AddAtomMolGraph(molecule, atom_idx, atom_type)
            for atom_idx in range(mol_graph.nb_atoms)
            if implicit_valences[atom_idx] > 0
            and not mol_graph.atom_charged_or_radical(atom_idx)
            for atom_type in Molecule.accepted_atoms
        ]
