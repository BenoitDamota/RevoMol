"""
Substitute an atom in the molecular graph.
"""

from copy import copy

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, Molecule, max_valence


class SubstituteAtomMolGraph(Action):
    """
    Substitute an atom in the molecular graph.
    In order to substitute an atom, the atom must be mutable.
    In addition, the substitution must be allowed in the
    accepted_substitutions dictionary.
    """

    accepted_substitutions: dict[str, list[str]]

    def __init__(self, molecule: Molecule, atom_idx: int, new_type: str) -> None:
        super().__init__(molecule)
        self.atom_idx: int = atom_idx
        self.new_type: str = new_type

    @classmethod
    def set_accepted_substitutions(
        cls, accepted_substitutions: dict[str, list[str]]
    ) -> None:
        """Set the accepted substitutions."""
        cls.accepted_substitutions = accepted_substitutions

    @classmethod
    def init_accepted_substitutions_from_accepted_atoms(cls) -> None:
        """Initialize the accepted substitutions from the accepted atoms."""
        cls.accepted_substitutions = {}
        for accepted_atom in Molecule.accepted_atoms:
            if accepted_atom != "H":
                curr_atom_subst: list[str] = []
                for other_accepted_atom in Molecule.accepted_atoms:
                    if other_accepted_atom != accepted_atom:
                        curr_atom_subst.append(other_accepted_atom)
                cls.accepted_substitutions[accepted_atom] = curr_atom_subst

    @override
    def apply(self) -> Molecule:
        mol_graph: MolecularGraph = self.molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph: MolecularGraph = copy(mol_graph)

        new_mol_graph.replace_atom(self.atom_idx, self.new_type)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Substitution atom caused an error.")
            raise e

        return Molecule(new_mol_graph.smiles)

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

        if not cls.accepted_substitutions:
            raise ValueError("Accepted substitutions not defined.")

        action_list: list[Action] = []

        max_valences = {
            atom_type: max_valence(atom_type)
            for atom_type in cls.accepted_substitutions
        }

        for atom in range(mol_graph.nb_atoms):
            atom_type = mol_graph.atom_type(atom)
            if not mol_graph.atom_mutability(atom):
                continue
            explicit_valence = mol_graph.atom_degree(atom, as_multigraph=True)
            for subst in cls.accepted_substitutions[atom_type]:
                if max_valences[subst] >= explicit_valence:
                    action_list.append(SubstituteAtomMolGraph(molecule, atom, subst))

        return action_list
