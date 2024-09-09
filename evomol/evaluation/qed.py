"""
Evaluation of population with QED score using RDKit implementation.

Bickerton, G. Richard, Gaia V. Paolini, Jérémy Besnard, Sorel Muresan, and
Andrew L. Hopkins.
Quantifying the Chemical Beauty of Drugs
Nature Chemistry 4, nᵒ 2 (février 2012): 90-98
https://doi.org/10.1038/nchem.1243

Based on RDKit implementation:
http://www.rdkit.org/new_docs/source/rdkit.Chem.QED.html
https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/QED.py

The choice of the weights for the QED score is the default one provided by RDKit
w=(0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)

To get the different scores :
# properties: Chem.QED.QEDproperties = Chem.QED.properties(
#     Chem.MolFromSmiles(mol_graph.canonical_smiles)
# )
# print(properties.ALERTS)
# print(properties.ALOGP)
# print(properties.AROM)
# print(properties.HBA)
# print(properties.HBD)
# print(properties.MW)
# print(properties.PSA)
# print(properties.ROTB)
"""

from rdkit.Chem import MolFromSmiles
from rdkit.Chem.QED import qed

from evomol.evaluation.evaluation import Function
from evomol.representation import MolecularGraph, Molecule


def qed_score(molecule: Molecule) -> float:
    """Calculate the quantitative estimate of drug-likeness (QED) score of a
    molecule.

    Args:
        molecule (Molecule): Molecule to evaluate

    Returns:
        float: QED score
    """
    mol_graph = molecule.get_representation(MolecularGraph)

    # QED is a score between 0 and 1
    qed_score_: float = qed(MolFromSmiles(mol_graph.canonical_smiles))

    return qed_score_


QED = Function("QED", qed_score)
