"""
This module contains the abstract class Evaluation.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Callable

from typing_extensions import override

from evomol.representation import Molecule


class EvaluationError(Exception):
    """Evaluation error."""

    def __init__(self, message: str) -> None:
        """Initialize the error.

        Args:
            message (str): The error message.
        """
        super().__init__(message)


class Evaluation(ABC):
    """
    Abstract class for evaluations.
    """

    def __init__(self, name: str) -> None:
        """Initialize the evaluation.

        Args:
            name (str): The name of the evaluation.
        """
        self._name: str = name

    @property
    def name(self) -> str:
        """Return the name of the evaluation.

        Returns:
            str: The name of the evaluation.
        """
        return self._name

    def evaluate(self, molecule: Molecule) -> Any:
        """Evaluate the property of the molecule. If the value is already
        computed, return it. Otherwise, compute it and store it in the molecule.

        Args:
            molecule (Molecule): Molecule to evaluate.

        Returns:
            object: The value of the evaluation.
        """
        value: Any = molecule.value(self.name)
        if value is not None:
            return value

        value = self._evaluate(molecule)

        molecule.set_value(self.name, value)
        return value

    @abstractmethod
    def _evaluate(self, molecule: Molecule) -> Any:
        """Evaluate the property of the molecule.

        Args:
            molecule (Molecule): Molecule to evaluate.

        Returns:
            object: The value of the evaluation.
        """


class Function(Evaluation):
    """Class to define an evaluation based on a function."""

    def __init__(self, name: str, function: Callable[[Molecule], Any]) -> None:
        super().__init__(name)
        self.function = function

    @override
    def _evaluate(self, molecule: Molecule) -> Any:
        return self.function(molecule)


def is_valid_molecule(molecule: Molecule, evaluations: list[Evaluation]) -> bool:
    """Check if a molecule is valid according to a list of evaluations.

    Args:
        molecule (Molecule): Molecule to check
        evaluations (list[Evaluation]): filter and evaluations to apply

    Returns:
        bool: True if the molecule is valid
    """
    for eval_ in evaluations:
        try:
            eval_.evaluate(molecule)
        except EvaluationError:
            return False
    return True
