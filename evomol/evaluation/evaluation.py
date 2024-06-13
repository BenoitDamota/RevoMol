from abc import ABC, abstractmethod
from typing import Any

from typing_extensions import override


from evomol.representation.molecule import Molecule


class Evaluation(ABC):

    def __init__(self, name: str):
        self._name = name

    @property
    def name(self):
        return self._name

    def evaluate(self, molecule: Molecule) -> dict[str, Any]:
        value = molecule.value(self.name)
        if value is not None:
            return value
        value = self._evaluate(molecule)
        molecule.set_value(self.name, value)
        return value

    @abstractmethod
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        pass
