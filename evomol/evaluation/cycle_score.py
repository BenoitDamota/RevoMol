"""
Cycle score evaluation.
"""

from typing import Any

import networkx as nx
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Molecule


class CycleScore(Evaluation):
    """
    Evaluation of molecule with cycle score.
    """

    def __init__(self) -> None:
        super().__init__("Cycle_score")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        mol_graph = molecule.get_representation(MolecularGraph)

        # normalization constants
        # statistics from 250k_rndm_zinc_drugs_clean.smi
        cycle_mean = -0.0485696876403053
        cycle_std = 0.2860212110245455

        # cycle score
        cycles = nx.cycle_basis(nx.Graph(mol_graph.adjacency_matrix))
        # or
        cycles = mol_graph.mol.GetRingInfo().AtomRings()

        cycle_length = max(len(cycle) for cycle in cycles) if cycles else 0

        cycle_length = max(cycle_length - 6, 0)

        cycle_score = -cycle_length

        normalized_cycle = (cycle_score - cycle_mean) / cycle_std

        return {
            "cycle_score": cycle_score,
            "normalized_cycle_score": normalized_cycle,
        }
