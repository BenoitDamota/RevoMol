"""
Cut an atom from the molecular graph.
"""

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, ActionSpace, Molecule


class CutAtomMolGraph(Action):
    """
    Cut an atom from the molecular graph.
    """

    def __init__(self, atom_to_cut: int, atom1_to_bond: int, atom2_to_bond: int):
        # index of the atom to cut
        self.atom_to_cut: int = atom_to_cut
        # index of the atom to bond to the new atom
        self.atom1_to_bond: int = atom1_to_bond
        self.atom2_to_bond: int = atom2_to_bond

    @override
    def apply(self, molecule: Molecule) -> Molecule:
        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None

        new_mol_graph: MolecularGraph = MolecularGraph(mol_graph.smiles)
        # remove bonds
        new_mol_graph.set_bond(self.atom_to_cut, self.atom1_to_bond, 0)
        new_mol_graph.set_bond(self.atom_to_cut, self.atom2_to_bond, 0)

        # create new bond
        new_mol_graph.set_bond(self.atom1_to_bond, self.atom2_to_bond, 1)

        new_mol_graph.remove_atom(self.atom_to_cut)

        try:
            new_mol_graph.update_representation()
        except Exception as e:
            print("Cut caused an error.")
            raise e

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return (
            "CutAtomMolGraph("
            f"{self.atom_to_cut}, {self.atom1_to_bond}, {self.atom2_to_bond})"
        )


class ActionSpaceCutAtomMolGraph(ActionSpace):
    """
    List possible actions on molecular graphs to cut an atom.
    """

    @override
    def list_actions(self, molecule: Molecule) -> list[Action]:
        """List possible actions to cut an atom from the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        action_list: list[Action] = []

        formal_charge_vector = mol_graph.formal_charge_vector()

        for atom_to_cut in range(mol_graph.nb_atoms):
            if not mol_graph.atom_mutability(atom_to_cut):
                continue
            bonds_to = [
                atom2_idx
                for atom2_idx in range(mol_graph.nb_atoms)
                if mol_graph.bond_type_num(atom_to_cut, atom2_idx) > 0
            ]
            if len(bonds_to) != 2:
                continue

            if (
                mol_graph.bond_type_num(bonds_to[0], bonds_to[1]) == 0
                and formal_charge_vector[bonds_to[0]] == 0
                and formal_charge_vector[bonds_to[1]] == 0
            ):
                action_list.append(
                    CutAtomMolGraph(atom_to_cut, bonds_to[0], bonds_to[1])
                )

        return action_list
