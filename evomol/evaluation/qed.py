"""
Evaluation of population with QED score using RDKit implementation.
"""

from typing import Any

from rdkit import Chem
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation import MolecularGraph, Molecule


class QED(Evaluation):
    """
    Evaluation of population with QED score using RDKit implementation.
    Bickerton, G. Richard, Gaia V. Paolini, Jérémy Besnard, Sorel Muresan, and
    Andrew L. Hopkins.
    Quantifying the Chemical Beauty of Drugs
    Nature Chemistry 4, nᵒ 2 (février 2012): 90-98
    https://doi.org/10.1038/nchem.1243
    """

    def __init__(self) -> None:
        super().__init__("QED")

    @override
    def _evaluate(self, molecule: Molecule) -> Any:

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        score = Chem.QED.qed(Chem.MolFromSmiles(mol_graph.canonical_smiles))

        return {
            "QED_score": score,
        }
