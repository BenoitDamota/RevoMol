"""
This script shows an example of a new action to add an atom to a molecule using
the SMILES representation.
The action is defined in the AddAtomSMILES class that inherits from the Action
class.
The action is applied to the molecule by adding an atom of a given type at a
position in the SMILES string.

This is not a complete and fully functional example, but a simple example to
show how to create a new action.
"""

import os
import sys

from typing_extensions import override

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# pylint: disable=wrong-import-position, import-error

from evomol.action import Action, ActionError, pprint_action_space_and_apply
from evomol.representation import SMILES, MolecularGraph, Molecule

# define the new action to add an atom to the molecule using SMILES


class AddAtomSMILES(Action):
    """Add an atom to the molecule using SMILES."""

    # the constructor of the class takes the inputs to apply the action

    def __init__(self, molecule: Molecule, position: int, atom_type: str) -> None:
        """Initialize the action.

        Args:
            molecule (Molecule): Molecule to which the atom is added
            position (int): position of the atom to add in the SMILES string
            atom_type (str): type of the new atom
        """
        super().__init__(molecule)

        # position of the atom in the SMILES string
        self.position: int = position
        # type of the new atom
        self.atom_type: str = atom_type

    # the __eq__, __hash__, __repr__, and representation_name are used
    # to compare actions or print

    @override
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AddAtomSMILES):
            return False
        return (
            self.molecule == other.molecule
            and self.position == other.position
            and self.atom_type == other.atom_type
        )

    @override
    def __hash__(self) -> int:
        return hash(self.__repr__())

    def __repr__(self) -> str:
        return f"AddAtomSMILES({self.molecule}, {self.position}, {self.atom_type})"

    @override
    def representation_name(self) -> str:
        return SMILES.class_name()

    # the _apply method is used to apply the action to the molecule,
    # it is called by the apply method of the Action class

    @override
    def _apply(self) -> Molecule:
        """Apply the action to the molecule.

        Raises:
            ActionError: Error if the action cannot be applied.

        Returns:
            Molecule: Molecule after the action.
        """
        # get the SMILES representation
        smiles: str = self.molecule.get_representation(SMILES).str_id

        # add the atom to the SMILES string
        new_smiles: str = (
            smiles[: self.position] + self.atom_type + smiles[self.position :]
        )

        # update the representation of the molecule using the MolecularGraph object
        try:
            new_mol_graph: MolecularGraph = MolecularGraph(new_smiles)
        except Exception as e:
            # raise an error if the new molecular graph cannot be converted to
            # a molecule
            raise ActionError(self, new_smiles, repr(e)) from e

        return Molecule(new_mol_graph.canonical_smiles)

    # the list_actions method is used to list all possible actions to add an atom
    # to the molecule

    @override
    @classmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List possible actions to add an atom to the SMILES.

        Three cases are considered:

        - if the molecule size is equal or higher than Molecule.max_heavy_atoms,
          no action is possible

        - if the molecule is empty, any atom can be added

        - otherwise, for each atom, if the next characters are letters from the
            accepted atoms, an action can be performed for each possible atom type

        Args:
            molecule (Molecule): Molecule to which the atom is added

        Returns:
            list[Action]: list of possible actions to add an atom to the SMILES
        """

        smiles: str = molecule.get_representation(SMILES).str_id

        accepted_letters = set(Molecule.accepted_atoms) | set(
            letter.lower() for letter in Molecule.accepted_atoms
        )

        # count the number of atoms in the molecule
        nb_atoms = 0
        for letter in smiles:
            if letter in accepted_letters:
                nb_atoms += 1

        # if the molecule is too big, we cannot add an atom
        if nb_atoms >= Molecule.max_heavy_atoms:
            return []

        # if the molecule is empty, we can add any atom
        if not smiles:
            return [
                AddAtomSMILES(molecule, 0, atom_type)
                for atom_type in Molecule.accepted_atoms
            ]

        # otherwise, for each atom, if the next characters are letters from the
        # accepted atoms, an action can be performed for each possible atom type
        return [
            AddAtomSMILES(molecule, position, atom_type)
            for position in range(len(smiles) + 1)
            for atom_type in Molecule.accepted_atoms
            if (position == 0 and smiles[position] in accepted_letters)
            or (position == len(smiles) and smiles[position - 1] in accepted_letters)
            or (
                smiles[position - 1] in accepted_letters
                and smiles[position] in accepted_letters
            )
        ]


def main() -> None:
    """Example of using the AddAtomSMILES action."""

    # initialize the parameters
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [SMILES]
    Molecule.max_heavy_atoms = 6

    Molecule.accepted_atoms = ["C", "O", "N", "S"]

    # initialize the action space for SMILES representation
    SMILES.action_space = [
        AddAtomSMILES,
    ]

    # create a molecule
    molecules = [
        Molecule(""),
        Molecule("C"),
        Molecule("CC"),
        Molecule("CCCCCC"),
    ]

    for molecule in molecules:
        # list possible actions
        possible_actions = molecule.compute_possible_actions()

        # print the possible actions
        pprint_action_space_and_apply(possible_actions, molecule)
        print()


if __name__ == "__main__":
    main()
