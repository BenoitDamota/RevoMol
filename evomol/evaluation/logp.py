from typing import Any

from typing_extensions import override
from rdkit.Chem import Descriptors

from evomol.evaluation import Evaluation
from evomol.representation.molecule import Molecule
from evomol.representation.molecular_graph import MolecularGraph


class LogP(Evaluation):

    def __init__(self):
        super().__init__("logP")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        # normalization constants
        # statistics from 250k_rndm_zinc_drugs_clean.smi
        logP_mean = 2.4570953396190123
        logP_std = 1.434324401111988

        mol_graph = molecule.get_representation(MolecularGraph)

        log_p = Descriptors.MolLogP(mol_graph.mol)

        normalized_log_p = (log_p - logP_mean) / logP_std

        return {
            "logP": log_p,
            "normalized_logP": normalized_log_p,
        }
