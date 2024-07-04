"""
List Morgan fingerprints that are not in the reference database.
ECFP4 are equivalent to Morgan fingerprints with radius 2.
Filter out molecules that contain unknown Morgan fingerprints.

Inspired on the work of Patrick Walters :
https://github.com/PatWalters/silly_walks

TODO add the reference to the paper
"""

import os

from rdkit import Chem
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation, EvaluationError
from evomol.representation import MolecularGraph, Molecule


class UnknownECFP(Evaluation):
    """
    Count Morgan fingerprints that are not in the reference database.
    ECFP4 are equivalent to Morgan fingerprints with radius 2.
    Filter out molecules that contain unknown Morgan fingerprints.
    """

    def __init__(
        self,
        path_db: str = os.path.join("external_data", "ecfp4_ChEMBL.txt"),
        radius: int = 2,
    ):
        """
        Init FilterECFP with the path to the reference database and the radius.

        Use file "external_data/ecfp4_ChEMBL_ZINC.txt" for ZINC database.

        Args:
            path_db (str): file that contains the reference list of ECFC4
                fingerprints keys
            radius (int, optional): radius if the ECFP fingerprint
                (2 for ECFP4, 1 for ECFP2). Defaults to 2.
        """
        super().__init__("UnknownECFP")

        self.radius = radius

        with open(path_db, encoding="utf8") as f:
            self.ecfp_list: frozenset[int] = frozenset(
                int(line.strip()) for line in f.readlines()
            )

    @override
    def _evaluate(self, molecule: Molecule) -> int:

        mol = Chem.MolFromSmiles(
            molecule.get_representation(MolecularGraph).canonical_smiles
        )

        if mol is None:
            return 1

        # get the Morgan fingerprints
        fingerprints = Chem.AllChem.GetMorganFingerprint(
            mol, radius=self.radius
        ).GetNonzeroElements()

        # count the number of ecfp that are not in the reference database
        nb_unknown_ecfp = 0
        for ecfp in fingerprints:
            if ecfp not in self.ecfp_list:
                nb_unknown_ecfp += 1

        return nb_unknown_ecfp


class FilterUnknownECFP(Evaluation):
    """Filter molecules with too many unknown ECFP."""

    def __init__(self, threshold: int = 0):
        super().__init__("FilterUnknownECFP")
        self.threshold = threshold

    @override
    def _evaluate(self, molecule: Molecule) -> bool:
        if molecule.value("UnknownECFP") > self.threshold:
            raise EvaluationError(
                "The molecule contains more unknown ECFP than allowed "
                f"({self.threshold})."
            )
        return True
