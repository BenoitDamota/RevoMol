import os
import sys

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


from evomol import evaluation as evaluator
from evomol.representation import SMILES, MolecularGraph, Molecule

Molecule.id_representation_class = SMILES
Molecule.representations_class = [MolecularGraph]

eval_name = "chembl"  # chembl chembl_zinc

evaluations = [
    evaluator.UnknownGCF(path_db="external_data/gcf1.txt", name="chembl"),
    evaluator.UnknownECFP(
        path_db="external_data/ecfp4_ChEMBL.txt", radius=2, name="chembl"
    ),
    evaluator.UnknownGCF(path_db="external_data/gcf2.txt", name="chembl_zinc"),
    evaluator.UnknownECFP(
        path_db="external_data/ecfp4_ChEMBL_ZINC.txt", radius=2, name="chembl_zinc"
    ),
]

smiles = "O=S(=O)(O)c1ccccc1c1ccsn1"
mol = Molecule(smiles)

print("canonical no aromatic SMILES:", mol)
print(
    "canonical aromatic SMILES   :",
    mol.get_representation(MolecularGraph).aromatic_canonical_smiles,
)

for eval_ in evaluations:
    print(eval_.name, eval_.evaluate(mol))
