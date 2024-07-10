"""
Action module.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import evomol.representation.molecule as molecule_module


class ActionContext(ABC):
    """
    Context of an action.

    The context is used to define the environment in which the action is
    applied.
    """

    def __init__(self, molecule: molecule_module.Molecule) -> None:
        self.molecule: molecule_module.Molecule = molecule

    @abstractmethod
    def get(self) -> object:
        """Return the context of the action."""


class Action(ABC):
    """
    Action to apply to a molecule representation.

    Apply function creates a new molecule based on the old one and apply the
    action to it.
    """

    def __init__(self, molecule: molecule_module.Molecule) -> None:
        self.molecule: molecule_module.Molecule = molecule

    def apply(self) -> molecule_module.Molecule:
        """Apply the action to the molecule."""
        # remove the action from the molecule possible actions
        self.molecule.add_applied_action(self)
        return self._apply()

    @abstractmethod
    def _apply(self) -> molecule_module.Molecule:
        """Apply the action"""

    @classmethod
    @abstractmethod
    def list_actions(cls, molecule: molecule_module.Molecule) -> list[Action]:
        """List all possible actions for the molecule"""

    @classmethod
    def class_name(cls) -> str:
        """Return the class name of the action."""
        return cls.__name__

    @abstractmethod
    def __eq__(self, value: object) -> bool:
        pass

    @abstractmethod
    def representation_name(self) -> str:
        """Return the name in str of the representation used in the action."""

    @abstractmethod
    def __hash__(self) -> int:
        """Return the hash of the action."""
        return hash(self.__repr__())


class ActionError(Exception):
    """Error raised when an action is not applicable to a molecule."""

    def __init__(self, action: Action, new_smiles: str, message: str) -> None:
        """Initialize the mutation error.

        Args:
            action (type[Action]): Action that caused the error
            new_smiles (str): SMILES after the mutation
            message (str): message of the error
        """
        self.message: str = (
            f"Action {action} caused an error. "
            f"New SMILES: {new_smiles}. "
            f"Message: {message}"
        )
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


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
