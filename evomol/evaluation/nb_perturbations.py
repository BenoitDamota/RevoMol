"""
This module contains the evaluation strategy that computes the number of
perturbations that were already applied to the molecule.
"""

from typing import Any

from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation import Molecule


class NPerturbations(Evaluation):
    """
    Strategy that computes the count of perturbations that were already applied
    to the solutions.
    If the parameters of EvoMol are set such as a mutation is a single
    perturbation on the molecular graph, then this value is equivalent to the
    number of previous mutations.
    """

    def __init__(self) -> None:
        super().__init__("NPerturbations")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        return {
            # "count": molecule.nb_perturbations,
        }
