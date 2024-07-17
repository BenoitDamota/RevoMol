"""
attempt to replace the tanimoto evaluation function from GuacaMol
but the results are not the same.
should explore the implementation in GuacaMol
"""

from guacamol.common_scoring_functions import TanimotoScoringFunction

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


target_smiles = "CCO"

scorer = TanimotoScoringFunction(
    target_smiles,
    fp_type="ECFP4",
)


canonical_red = Chem.CanonSmiles(target_smiles)
mol_red = Chem.MolFromSmiles(canonical_red)
fp_red = FingerprintMols.FingerprintMol(mol_red)
fp_morgan_red = Chem.AllChem.GetMorganFingerprint(mol_red, 2)

for smiles in ["NNNNNN", "CCOCC", "CCO"]:
    canonical = Chem.CanonSmiles(smiles)
    mol = Chem.MolFromSmiles(canonical)
    fp = FingerprintMols.FingerprintMol(mol)
    fp_morgan = Chem.AllChem.GetMorganFingerprint(mol, 2)

    score_data = DataStructs.FingerprintSimilarity(fp, fp_red)

    score_data_morgan = DataStructs.DiceSimilarity(fp_morgan, fp_morgan_red)

    score_guaca = scorer.score(canonical)

    print(
        f"Score for {smiles:10}: {score_data:.2f} (data) vs {score_data_morgan:.2f} (morgan) vs {score_guaca:.2f} (guacamol)"
    )
