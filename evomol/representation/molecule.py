"""
Representation of a molecule and its operations.

The following page can be useful :
    https://www.rdkit.org/docs/cppapi/classRDKit_1_1RWMol.html
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar, TypeVar, Type

from rdkit import Chem

# Maximum valence for sulfur
SULFUR_MAX_VALENCE = 6


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


def pprint_action_space(
    action_list: dict[str, dict[str, list[Action]]],
) -> None:
    """Pretty print the action space of a molecule."""
    nb_actions = 0
    print(f"{len(action_list)} representations: ")
    for representation, action_spaces in action_list.items():
        print(repr(representation))
        print(f"\t{len(action_spaces.keys())} actions space :")
        for action_space, actions in action_spaces.items():
            print(f"\t{action_space}")
            print(f"\t\t{len(actions)} actions :")
            for action in actions:
                print(f"\t\t{action}")
                nb_actions += 1
    print(f"Nb actions: {nb_actions}")


class ActionContext(ABC):
    """
    Context of an action.

    The context is used to define the environment in which the action is
    applied.
    """

    def __init__(self, molecule: Molecule) -> None:
        self.molecule: Molecule = molecule

    @abstractmethod
    def get(self) -> object:
        """Return the context of the action."""


class Action(ABC):
    """
    Action to apply to a molecule representation.

    Apply function creates a new molecule based on the old one and apply the
    action to it.
    """

    def __init__(self, molecule: Molecule) -> None:
        self.molecule: Molecule = molecule

    @abstractmethod
    def apply(self) -> Molecule:
        """Apply the action"""

    @classmethod
    @abstractmethod
    def list_actions(cls, molecule: Molecule) -> list[Action]:
        """List all possible actions for the molecule"""

    @classmethod
    def class_name(cls) -> str:
        """Return the class name of the action."""
        return cls.__name__


# TypeAction = TypeVar("TypeAction", bound=Action)


class MoleculeRepresentation(ABC):
    """
    Representation of a molecule.

    The purpose of this class is to represent a molecule in a way that can be
    used to generate neighbors (applied actions) or identify it.

    For each of the representations, the possible actions spaces are defined
    by the action_space attribute.
    This attribute must be defined in the child class to give the ability to
    list all possible actions for a molecule.
    To set the possible neighbors :
        YourMoleculeRepresentationClass.set_action_spaces(
            [
                Action1,
                Action2,
                ...,
            ]
        )

        # or
        YourMoleculeRepresentationClass.action_space = [
            Action1,
            Action2,
            ...,
        ]

    """

    # possible action space for the molecule representation
    action_space: ClassVar[list[Type[Action]]]

    def __init__(self, str_id: str):
        """init the molecule representation from a string."""
        self.str_id: str = str_id

    @classmethod
    def set_action_spaces(cls, action_space: list[Type[Action]]) -> None:
        """Set the possible action space for the molecule representation."""
        cls.action_space = action_space

    def list_all_actions(self, molecule: Molecule) -> dict[str, list[Action]]:
        """
        List all possible Action for each Action of the molecule representation
        """
        return {
            action.class_name(): action.list_actions(molecule)
            for action in self.action_space
        }

    @abstractmethod
    def representation(self) -> str:
        """String used to identify the molecule with this representation"""
        raise NotImplementedError

    def __repr__(self) -> str:
        return self.__name__

    @classmethod
    def class_name(self) -> str:
        """Return the class name of the representation."""
        return self.__name__

    def __str__(self) -> str:
        return self.representation()


TypeMoleculeRepresentation = TypeVar(
    "TypeMoleculeRepresentation", bound=MoleculeRepresentation
)


class Molecule:
    """
    Molecule with its various representations.

    The id_representation is used to define the molecule in a unique way for
    the whole system (id for cache, logs, ...).
    The representations are used to generate neighbors or to compute properties
    of the molecule.
    """

    id_representation_class: type[MoleculeRepresentation]
    representations_class: list[type[MoleculeRepresentation]]
    max_heavy_atoms: int
    accepted_atoms: list[str] = ["C", "O", "N", "F"]

    def __init__(self, str_id: str):
        """init the molecule from a string and create its various representations.

        Args:
            str_id (str): string used to identify the molecule(SMILES, InChI, ...)
        """
        self.id_representation: MoleculeRepresentation = (
            Molecule.id_representation_class(str_id)
        )
        self.representations: list[MoleculeRepresentation] = [
            representation_class(str_id)
            for representation_class in Molecule.representations_class
        ]

    def list_available_actions_space(self) -> dict[str, list[Type[Action]]]:
        """List all available actions space for each representation."""
        return {
            representation.class_name(): representation.action_space
            for representation in self.representations
            if representation.action_space
        }

    def list_all_possible_actions(self) -> dict[str, dict[str, list[Action]]]:
        """List all possible actions for the molecule."""
        # make sure that the possible action is not empty
        possible_actions = {}
        for representation in self.representations:
            current_actions = {}
            for action in representation.action_space:
                list_actions = action.list_actions(self)
                if list_actions:
                    current_actions[action.class_name()] = list_actions
            if current_actions:
                possible_actions[representation.class_name()] = current_actions
        return possible_actions

        # return {
        #     representation: {
        #         action_space: action_space.list_actions(self)
        #         for action_space in representation.list_actions_space(self)
        #     }
        #     for representation in self.representations
        # }

    def get_representation(
        self, representation_class: Type[TypeMoleculeRepresentation]
    ) -> TypeMoleculeRepresentation:
        """Return the representation of the molecule for the given class."""
        for representation in self.representations:
            if type(representation) is representation_class:
                return representation
        raise ValueError(f"No representation found for {representation_class}")

    @classmethod
    def get_action_space(
        cls,
        representation_class: Type[MoleculeRepresentation],
        action_space_class: Type[Action],
    ) -> Type[Action]:
        """Return the action space for the given representation class."""
        for representation in cls.representations_class:
            if representation == representation_class:
                for action_space in representation.action_space:
                    if action_space is action_space_class:
                        return action_space
        raise ValueError(
            f"No action space found for {action_space_class} in {representation_class}"
        )

    def __repr__(self) -> str:
        return repr(self.id_representation)

    def __str__(self) -> str:
        return str(self.id_representation)


# voir pour faire des tests avec des neighborhood strategies et appliquer k actions
