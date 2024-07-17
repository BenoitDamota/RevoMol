"""
Rediscovery score based on the implementation of GuacaMol.

Removed to avoid guacamol dependency.
Look at scripts/tanimoto.py for an attempt to use only RDKit to calculate the
similarity.

Nathan Brown et al.
GuacaMol: Benchmarking Models for de Novo Molecular Design
Journal of Chemical Information and Modeling 59, no. 3 (March 25, 2019): 1096–1108
https://doi.org/10.1021/acs.jcim.8b00839
"""

# from guacamol.common_scoring_functions import TanimotoScoringFunction
# from typing_extensions import override

# from evomol.evaluation.evaluation import Evaluation
# from evomol.representation import MolecularGraph, Molecule


# class Rediscovery(Evaluation):
#     """
#     Rediscovery score based on the implementation of GuacaMol
#     """

#     def __init__(self, target_smiles: str) -> None:
#         super().__init__("rediscovery_" + target_smiles)
#         self.scorer = TanimotoScoringFunction(
#             target_smiles,
#             fp_type="ECFP4",
#         )

#     @override
#     def _evaluate(self, molecule: Molecule) -> float:
#         mol_graph = molecule.get_representation(MolecularGraph)

#         return self.scorer.score(mol_graph.aromatic_canonical_smiles)
