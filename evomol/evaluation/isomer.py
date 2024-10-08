"""
Isomer score based on the implementation of GuacaMol.

Removed to avoid guacamol dependency.

Nathan Brown et al.
GuacaMol: Benchmarking Models for de Novo Molecular Design
Journal of Chemical Information and Modeling 59, no. 3
(March 25, 2019): 1096–1108
https://doi.org/10.1021/acs.jcim.8b00839
"""

# from guacamol.common_scoring_functions import IsomerScoringFunction
# from typing_extensions import override

# from evomol.evaluation.evaluation import Evaluation
# from evomol.representation import MolecularGraph, Molecule


# class Isomer(Evaluation):
#     """
#     Isomer score based on the implementation of GuacaMol
#     """

#     def __init__(self, formula: str):
#         super().__init__("isomer_guacamol")
#         self.scorer = IsomerScoringFunction(formula)

#     @override
#     def _evaluate(self, molecule: Molecule) -> float:
#         mol_graph = molecule.get_representation(MolecularGraph)

#         return self.scorer.score(mol_graph.aromatic_canonical_smiles)
