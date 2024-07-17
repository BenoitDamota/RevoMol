"""
Rediscovery score based on the implementation of GuacaMol.

Removed to avoid guacamol dependency.

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


# This code is an attempt to use only RDKit to calculate the similarity but
# the results are not the same as the GuacaMol implementation.
# should explore the implementation in GuacaMol

# from rdkit import Chem
# from rdkit import DataStructs
# from rdkit.Chem.Fingerprints import FingerprintMols

# smiles_to_rediscover = "CCO"
# canonical_red = Chem.CanonSmiles(smiles_to_rediscover)
# mol_red = Chem.MolFromSmiles(canonical_red)
# fp_red = FingerprintMols.FingerprintMol(mol_red)


# smiles = "NNNNNN"
# smiles = "CCOCC"
# smiles = "CCO"
# canonical = Chem.CanonSmiles(smiles)
# mol = Chem.MolFromSmiles(canonical)
# fp = FingerprintMols.FingerprintMol(mol)

# score = DataStructs.FingerprintSimilarity(fp, fp_red)
# score
