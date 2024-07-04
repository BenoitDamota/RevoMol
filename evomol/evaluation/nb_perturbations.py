"""
This module contains the evaluation strategy that computes the number of
perturbations that were already applied to the molecule.
"""

from evomol.evaluation.evaluation import Function
from evomol.representation import Molecule


def n_perturbations(molecule: Molecule) -> int:
    """Calculate the number of perturbations that were already applied to the
    molecule.

    Args:
        molecule (Molecule): Molecule to evaluate

    Returns:
        int: Number of perturbations
    """
    _ = molecule
    return 0  # TODO molecule.nb_perturbations


NPerturbations = Function("NPerturbations", n_perturbations)
