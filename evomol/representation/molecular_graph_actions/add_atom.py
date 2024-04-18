"""
Add an atom to a molecular graph.
"""

from typing_extensions import override

from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Action, ActionSpace, Molecule


class AddAtomMolGraph(Action):
    """
    Add an atom to the molecular graph.
    """

    def __init__(self, index_atom: int, atom_type: str) -> None:
        # index of the atom to bond to the new atom
        self.index_atom: int = index_atom
        # type of the new atom
        self.atom_type: str = atom_type

    @override
    def apply(self, molecule: Molecule) -> Molecule:
        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        new_mol_graph = MolecularGraph(mol_graph.smiles)

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

        return Molecule(new_mol_graph.canonical_smiles)

    def __repr__(self) -> str:
        return f"AddAtomMolGraph({self.index_atom}, {self.atom_type})"


class ActionSpaceAddAtomMolGraph(ActionSpace):
    """
    List possible actions on molecular graphs to add an atom.
    """

    def __init__(self, accepted_atoms: list[str], allow_bonding: bool = False):
        """Initialize the Add Atom neighborhood operation.


        Args:
            allow_bonding (bool, optional):
                If True, atoms can be added either without bonding or with
                bonding to an atom that has free electrons and is connectable.
                Defaults to False.
        """
        self.accepted_atoms: list[str] = accepted_atoms
        self.allow_bonding: bool = allow_bonding

    @override
    def list_actions(self, molecule: Molecule) -> list[Action]:
        """List possible actions to add an atom to the molecular graph."""

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)
        assert mol_graph is not None
        # if not self.allow_bonding: TODO remove
        #     if mol_graph.nb_atoms >= Molecule.max_heavy_atoms:
        #         return []
        #     # all atoms are valid if the molecule is not at its maximum size
        #     return [AddAtomMolGraph(0, atom_type)
        # for atom_type in self.accepted_atoms]

        action_list: list[Action] = []

        if mol_graph.nb_atoms == 0:
            # possibility of adding an unconnected atom if the molecule is empty
            return [AddAtomMolGraph(0, atom_type) for atom_type in self.accepted_atoms]

        # possibility of adding an atom with a bond
        if mol_graph.nb_atoms < Molecule.max_heavy_atoms:
            free_electrons_vector = mol_graph.free_electrons_vector()
            formal_charge_vector = mol_graph.formal_charge_vector()

            # the new atom can be bonded to an existing atom if it has free
            # electrons and has no formal charge
            for atom_idx in range(mol_graph.nb_atoms):
                action_list.extend(
                    AddAtomMolGraph(atom_idx, atom_type)
                    for atom_type in self.accepted_atoms
                    if free_electrons_vector[atom_idx] > 0
                    and formal_charge_vector[atom_idx] == 0
                )

        return action_list
