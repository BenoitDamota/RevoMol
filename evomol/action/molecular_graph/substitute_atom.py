"""
Substitute an atom in the molecular graph.
"""

from rdkit import Chem
from typing_extensions import override

from evomol.action import Action
from evomol.representation import SULFUR_MAX_VALENCE, MolecularGraph, Molecule

from .action_molecular_graph import ActionMolGraph


def max_valence(atom: str) -> int:
    """Return the maximum valence of an atom.
    Valence for sulfur is set with the SULFUR_MAX_VALENCE constant
    (set to 6 by default).

    Args:
        atom (str): atom symbol

    Returns:
        int: valence of the atom
    """
    if atom == "S":
        return SULFUR_MAX_VALENCE

    valence: int = Chem.GetPeriodicTable().GetDefaultValence(
        Chem.GetPeriodicTable().GetAtomicNumber(atom)
    )
    return valence


class SubstituteAtomMG(ActionMolGraph):
    """
    Substitute an atom in the molecular graph.
    """

    def __init__(self, molecule: Molecule, atom_idx: int, new_type: str) -> None:
        super().__init__(molecule)
        self.atom_idx: int = atom_idx
        self.new_type: str = new_type

    @override
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SubstituteAtomMG):
            return False
        return (
            self.molecule == other.molecule
            and self.atom_idx == other.atom_idx
            and self.new_type == other.new_type
        )

    @override
    def __hash__(self) -> int:
        return hash(self.__repr__())

    @override
    def apply_action(self, new_mol_graph: MolecularGraph) -> None:
        new_mol_graph.replace_atom(self.atom_idx, self.new_type)

    def __repr__(self) -> str:
        return (
            f"SubstituteAtomMolGraph({self.molecule}, "
            f"{self.atom_idx}, {self.new_type})"
        )

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to substitute an atom in the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        max_valences = {
            atom_type: max_valence(atom_type) for atom_type in Molecule.accepted_atoms
        }

        for atom in range(mol_graph.nb_atoms):
            atom_type = mol_graph.atom_type(atom)
            if mol_graph.atom_charged_or_radical(atom):
                continue
            explicit_valence = mol_graph.atom_degree(atom, as_multigraph=True)
            for subst in Molecule.accepted_atoms:
                if subst == atom_type:
                    continue
                if max_valences[subst] >= explicit_valence:
                    action_list.append(SubstituteAtomMG(molecule, atom, subst))

        return action_list
