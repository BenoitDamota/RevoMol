"""
This module contains the evaluation class for the logP property.
"""

from typing import Any

from rdkit.Chem import Descriptors
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Molecule


class LogP(Evaluation):
    """
    Evaluation class for the logP property.
    """

    def __init__(self) -> None:
        super().__init__("logP")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        # normalization constants
        # statistics from 250k_rndm_zinc_drugs_clean.smi
        log_p_mean = 2.4570953396190123
        log_p_std = 1.434324401111988

        mol_graph = molecule.get_representation(MolecularGraph)

        log_p = Descriptors.MolLogP(mol_graph.mol)

        normalized_log_p = (log_p - log_p_mean) / log_p_std

        return {
            "logP": log_p,
            "normalized_logP": normalized_log_p,
        }
