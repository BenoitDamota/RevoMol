"""
Cycle score evaluation.
"""

import networkx as nx

from evomol.evaluation.evaluation import Function
from evomol.representation import MolecularGraph, Molecule


def cycle_score(molecule: Molecule) -> float:
    """Compute the cycle score of a molecule.

    Args:
        molecule (Molecule): Molecule to evaluate.

    Returns:
        float: The cycle score of the molecule.
    """
    mol_graph = molecule.get_representation(MolecularGraph)

    # find cycles in the molecule
    cycles = nx.cycle_basis(nx.Graph(mol_graph.adjacency_matrix))

    # find maximum cycle length
    cycle_length = max(len(cycle) for cycle in cycles) if cycles else 0

    # only consider cycles with length greater than 6
    cycle_length = max(cycle_length - 6, 0)

    cycle_score_ = -cycle_length

    return cycle_score_


def zinc_normalized_cycle_score(molecule: Molecule) -> float:
    """Normalize the cycle score of a molecule according to the ZINC database.

    Args:
        molecule (Molecule): Molecule to evaluate.

    Returns:
        float: The normalized cycle score of the molecule.
    """
    cycle_score_: float = molecule.value("CycleScore")

    # normalization constants
    # statistics from 250k_rndm_zinc_drugs_clean.smi
    cycle_mean: float = -0.0485696876403053
    cycle_std: float = 0.2860212110245455
    normalized_cycle: float = (cycle_score_ - cycle_mean) / cycle_std

    return normalized_cycle


CycleScore = Function("CycleScore", cycle_score)

NormalizedCycleScore = Function(
    "zinc_normalized_cycle_score", zinc_normalized_cycle_score
)
