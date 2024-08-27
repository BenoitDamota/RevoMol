"""
SAscore evaluation of a molecule using RDKit implementation.


Ertl, Peter, et Ansgar Schuffenhauer.
Estimation of synthetic accessibility score of drug-like
molecules based on molecular complexity and fragment contributions
Journal of Cheminformatics 1, no 1 (10 juin 2009): 8.
https://doi.org/10.1186/1758-2946-1-8.
"""

from rdkit.Contrib.SA_Score import sascorer

from evomol.evaluation.evaluation import Function
from evomol.representation import MolecularGraph, Molecule


def sa_score(molecule: Molecule) -> float:
    """Calculate the synthetic accessibility score of a molecule.

    Args:
        molecule (Molecule): Molecule to evaluate

    Returns:
        float: Synthetic accessibility score
    """
    mol_graph = molecule.get_representation(MolecularGraph)

    # SAscore is a positive score between 1 and 10
    # 10 being the most complex molecule
    sa_score_value: float = sascorer.calculateScore(mol_graph.mol)

    return sa_score_value


def normalized_sa_score(molecule: Molecule) -> float:
    """Normalize the synthetic accessibility score of a molecule.

    Args:
        molecule (Molecule): Molecule to evaluate

    Returns:
        float: Normalized synthetic accessibility score
    """
    sa_score_value: float = molecule.value("SAscore")

    # normalization with 0 being the most complex molecule
    normalized_sa_score_value: float = 1 - (sa_score_value - 1) / 9

    return normalized_sa_score_value


def zinc_normalized_sa_score(molecule: Molecule) -> float:
    """Normalize the synthetic accessibility score of a molecule using the
    statistics from 250k_rndm_zinc_drugs_clean.smi.

    Args:
        molecule (Molecule): Molecule to evaluate

    Returns:
        float: Normalized synthetic accessibility score using ZINC statistics
    """
    sa_score_value: float = molecule.value("SAscore")

    # normalization constants
    # statistics from 250k_rndm_zinc_drugs_clean.smi
    sa_mean: float = -3.0525811293166134
    sa_std: float = 0.8335207024513095
    zinc_normalized_sa_score_value: float = (-sa_score_value - sa_mean) / sa_std

    return zinc_normalized_sa_score_value


SAScore = Function("SAScore", sa_score)


NormalizedSAScore = Function("NormalizedSAScore", normalized_sa_score)


ZincNormalizedSAScore = Function(
    "zinc_normalized_sa_score",
    zinc_normalized_sa_score,
)
