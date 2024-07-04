"""
This module contains the evaluation class for the logP property.
"""

from rdkit.Chem import Descriptors

from evomol.evaluation.evaluation import Function
from evomol.representation import MolecularGraph, Molecule


def log_p(molecule: Molecule) -> float:
    """Calculate the logP of a molecule.

    Args:
        molecule (Molecule): Molecule to evaluate

    Returns:
        float: logP value
    """
    mol_graph = molecule.get_representation(MolecularGraph)

    log_p_value: float = Descriptors.MolLogP(mol_graph.mol)

    return log_p_value


def zinc_normalized_log_p(molecule: Molecule) -> float:
    """Normalize the logP of a molecule using the statistics from
    250k_rndm_zinc_drugs_clean.smi.

    Args:
        molecule (Molecule): Molecule to evaluate

    Returns:
        float: Normalized logP using ZINC statistics
    """
    log_p_value: float = molecule.value("logP")

    # normalization constants
    # statistics from 250k_rndm_zinc_drugs_clean.smi
    log_p_mean: float = 2.4570953396190123
    log_p_std: float = 1.434324401111988

    normalized_log_p: float = (log_p_value - log_p_mean) / log_p_std

    return normalized_log_p


LogP = Function("logP", log_p)
ZincNormalizedLogP = Function("zinc_normalized_logP", zinc_normalized_log_p)
