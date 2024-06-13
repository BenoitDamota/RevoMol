from typing import Any

from typing_extensions import override
from guacamol.common_scoring_functions import TanimotoScoringFunction

from evomol.evaluation import Evaluation
from evomol.representation.molecule import Molecule
from evomol.representation.molecular_graph import MolecularGraph


class Rediscovery(Evaluation):
    """
    Rediscovery score based on the implementation of GuacaMol
    Nathan Brown et al.
    GuacaMol: Benchmarking Models for de Novo Molecular Design
    Journal of Chemical Information and Modeling 59, no. 3 (March 25, 2019): 1096–1108
    https://doi.org/10.1021/acs.jcim.8b00839
    """

    def __init__(self, target_smiles: str):
        super().__init__("rediscovery_" + target_smiles)
        self.scorer = TanimotoScoringFunction(
            target_smiles,
            fp_type="ECFP4",
        )

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        mol_graph = molecule.get_representation(MolecularGraph)

        return {
            "rediscovery_score": self.scorer.score(mol_graph.canonical_smiles),
        }
