"""
Substitute an atom in the molecular graph.
"""

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, ActionSpace, Molecule, max_valence


class SubstituteAtomMolGraph(Action):
    """
    Substitute an atom in the molecular graph.
    In order to substitute an atom, the atom must be mutable.
    In addition, the substitution must be allowed in the
    accepted_substitutions dictionary.
    """

    def __init__(self, atom_idx: int, new_type: str) -> None:
        self.atom_idx: int = atom_idx
        self.new_type: str = new_type

    @override
    def apply(self, molecule: Molecule) -> Molecule:
        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph = MolecularGraph(mol_graph.smiles)

        new_mol_graph.replace_atom(self.atom_idx, self.new_type)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Substitution atom caused an error.")
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return f"SubstituteAtomMolGraph({self.atom_idx}, {self.new_type})"


class ActionSpaceSubstituteAtomMolGraph(ActionSpace):
    """
    List possible actions on molecular graphs to substitute an atom.
    """

    def __init__(self, accepted_substitutions: dict[str, list[str]]) -> None:
        self.accepted_substitutions: dict[str, list[str]] = accepted_substitutions

    @override
    def list_actions(self, molecule: Molecule) -> list[Action]:
        """List possible actions to substitute an atom in the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        max_valences = {
            atom_type: max_valence(atom_type)
            for atom_type in self.accepted_substitutions
        }

        for atom in range(mol_graph.nb_atoms):
            atom_type = mol_graph.atom_type(atom)
            explicit_valence = mol_graph.atom_degree(atom, as_multigraph=True)
            for subst in self.accepted_substitutions[atom_type]:
                if max_valences[subst] >= explicit_valence:
                    action_list.append(SubstituteAtomMolGraph(atom, subst))

        return action_list
