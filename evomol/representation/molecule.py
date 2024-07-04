"""
Representation of a molecule and its operations.

The following page can be useful :
    https://www.rdkit.org/docs/cppapi/classRDKit_1_1RWMol.html
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, ClassVar, TypeVar

import evomol.action.action as action_module  # pylint: disable=cyclic-import

# Maximum valence for sulfur
SULFUR_MAX_VALENCE = 6


class MoleculeRepresentation(ABC):
    """
    Representation of a molecule.

    The purpose of this class is to represent a molecule in a way that can be
    used to generate neighbors (applied actions) or identify it.

    For each of the representations, the possible actions spaces are defined
    by the action_space attribute.
    This attribute must be defined in the child class to give the ability to
    list all possible actions for a molecule.
    """

    # possible action space for the molecule representation
    action_space: ClassVar[list[type[action_module.Action]]] = []

    def __init__(self, str_id: str):
        """init the molecule representation from a string."""
        self.str_id: str = str_id

    @classmethod
    def set_action_space(
        cls,
        action_space: list[type[action_module.Action]],
    ) -> None:
        """Set the possible action space for the molecule representation."""
        cls.action_space = action_space

    @classmethod
    def list_action_spaces(cls) -> list[type[action_module.Action]]:
        """List all possible actions for the molecule representation."""
        return cls.action_space

    @abstractmethod
    def representation(self) -> str:
        """String used to identify the molecule with this representation"""
        raise NotImplementedError

    def __repr__(self) -> str:
        return MoleculeRepresentation.__name__

    @classmethod
    def class_name(cls) -> str:
        """Return the class name of the representation."""
        return cls.__name__

    def __str__(self) -> str:
        return self.representation()

    @abstractmethod
    def __eq__(self, value: object) -> bool:
        raise NotImplementedError


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
    max_heavy_atoms: int = 38
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
        self._values: dict[str, object] = {}
        self.possible_actions: dict[str, dict[str, list[action_module.Action]]] = {}

    def list_available_actions_space(
        self,
    ) -> dict[str, list[type[action_module.Action]]]:
        """List all available actions space for each representation."""
        return {
            representation.class_name(): representation.action_space
            for representation in self.representations
            if representation.action_space
        }

    def list_all_possible_actions(
        self,
    ) -> dict[str, dict[str, list[action_module.Action]]]:
        """List all possible actions for the molecule."""
        # if the possible actions have already been computed
        if self.possible_actions != {} and self.possible_actions is not None:
            return self.possible_actions

        # if all the possible actions have been applied (all removed from the list)
        if self.possible_actions is None:
            return None

        # compute the possible actions
        for representation in self.representations:
            current_actions = {}
            for action in representation.action_space:
                list_actions = action.list_actions(self)
                if list_actions:
                    current_actions[action.class_name()] = list_actions
            if current_actions:
                self.possible_actions[representation.class_name()] = current_actions
        return self.possible_actions

        # return {
        #     representation: {
        #         action_space: action_space.list_actions(self)
        #         for action_space in representation.list_actions_space(self)
        #     }
        #     for representation in self.representations
        # }

    def nb_possible_actions(self) -> int:
        """Return the number of possible actions for the molecule."""
        if self.possible_actions is None:
            return 0
        if not self.possible_actions:
            self.list_all_possible_actions()

        count = 0
        for representation in self.possible_actions.values():
            for actions in representation.values():
                count += len(actions)
        return count

    def remove_action(self, action: action_module.Action) -> None:
        """Remove the action from the possible actions."""
        repr_name = ""
        for representation_name, actions_dict in self.possible_actions.items():
            if repr_name:
                break
            for action_name in actions_dict:
                if action.class_name() == action_name:
                    repr_name = representation_name
                    break
        if not repr_name:
            raise ValueError(f"Action {action} not found")
        self.possible_actions[repr_name][action.class_name()].remove(action)
        if not self.possible_actions[repr_name][action.class_name()]:
            self.possible_actions[repr_name].pop(action.class_name())
        if not self.possible_actions[repr_name]:
            self.possible_actions.pop(repr_name)
        if not self.possible_actions:
            self.possible_actions = None

    def get_representation(
        self, representation_class: type[TypeMoleculeRepresentation]
    ) -> TypeMoleculeRepresentation:
        """Return the representation of the molecule for the given class."""
        for representation in self.representations:
            if isinstance(representation, representation_class):
                return representation
        raise ValueError(f"No representation found for {representation_class}")

    def value(self, name: str) -> Any:
        """Return the value of the molecule."""
        return self._values.get(name, None)

    def set_value(self, name: str, value: object) -> None:
        """Set the value of the molecule."""
        self._values[name] = value

    @classmethod
    def list_representation_classes(cls) -> list[type[MoleculeRepresentation]]:
        """List all possible representations for the molecule."""
        return cls.representations_class

    def __repr__(self) -> str:
        return repr(self.id_representation)

    def __str__(self) -> str:
        return str(self.id_representation)

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, Molecule)
            and self.id_representation == value.id_representation
            and self.representations == value.representations
        )

    def __hash__(self) -> int:
        # TODO to correct, MolecularGraph may not always be the first representation
        return hash(self.representations[0].canonical_smiles)
