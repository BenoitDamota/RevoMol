"""
This module contains the evaluation class for the penalized logP property.
"""

from evomol.evaluation.evaluation import Function
from evomol.representation import Molecule


def penalized_log_p(molecule: Molecule) -> float:
    """Calculate the penalized logP of a molecule using the zinc statistics.

    Args:
        molecule (Molecule): Molecule to evaluate

    Returns:
        float: PLogP value
    """

    n_log_p: float = molecule.value("zinc_normalized_logP")
    n_sa_score: float = molecule.value("zinc_normalized_sa_score")
    n_cycle_score: float = molecule.value("zinc_normalized_cycle_score")

    return n_log_p + n_sa_score + n_cycle_score


PLogP = Function("PLogP", penalized_log_p)
