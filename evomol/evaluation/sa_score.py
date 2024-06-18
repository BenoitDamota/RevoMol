"""
Evaluation of molecule with SA score using RDKit implementation.
"""

from typing import Any

from rdkit.Contrib.SA_Score import sascorer
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Molecule


class SAScore(Evaluation):
    """
    Evaluation of molecule with SA score using RDKit implementation.
    """

    def __init__(self) -> None:
        super().__init__("SA_score")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        mol_graph = molecule.get_representation(MolecularGraph)

        sa_score = -sascorer.calculateScore(mol_graph.mol)

        # normalization constants
        # statistics from 250k_rndm_zinc_drugs_clean.smi
        sa_mean = -3.0525811293166134
        sa_std = 0.8335207024513095

        normalized_sa_score = (sa_score - sa_mean) / sa_std

        # Ertl, Peter, et Ansgar Schuffenhauer.
        # Estimation of synthetic accessibility score of drug-like
        # molecules based on molecular complexity and fragment contributions
        # Journal of Cheminformatics 1, no 1 (10 juin 2009): 8.
        # https://doi.org/10.1186/1758-2946-1-8.
        other_sa_score = 1 - (sa_score - 1) / 9

        return {
            "SA_score": sa_score,
            "normalized_SA_score": normalized_sa_score,
            "other_SA_score": other_sa_score,
        }
