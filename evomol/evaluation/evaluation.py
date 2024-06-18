"""
This module contains the abstract class Evaluation.
"""

from abc import ABC, abstractmethod

from evomol.representation import Molecule


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

    def evaluate(self, molecule: Molecule) -> object:
        """Evaluate the property of the molecule. If the value is already
        computed, return it. Otherwise, compute it and store it in the molecule.

        Args:
            molecule (Molecule): Molecule to evaluate.

        Returns:
            object: The value of the evaluation.
        """
        value: object = molecule.value(self.name)
        if value is not None:
            return value
        value = self._evaluate(molecule)
        molecule.set_value(self.name, value)
        return value

    @abstractmethod
    def _evaluate(self, molecule: Molecule) -> object:
        """Evaluate the property of the molecule.

        Args:
            molecule (Molecule): Molecule to evaluate.

        Returns:
            object: The value of the evaluation.
        """
