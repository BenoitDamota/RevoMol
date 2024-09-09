"""
List Morgan fingerprints that are not in the reference database.
ECFP4 are equivalent to Morgan fingerprints with radius 2.
Filter out molecules that contain unknown Morgan fingerprints.

Inspired on the work of Patrick Walters :
https://github.com/PatWalters/silly_walks
"""

import os

from rdkit import Chem
from rdkit.Chem import AllChem
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation, EvaluationError
from evomol.representation import MolecularGraph, Molecule


def list_ecfp(molecule: Molecule, radius: int = 2) -> list[int]:
    """
    List the ECFP fingerprints of a molecule.

    Args:
        molecule (Molecule): molecule to evaluate
        radius (int, optional): radius of the ECFP fingerprint. Defaults to 2.

    Returns:
        list[int]: list of ECFP fingerprints
    """
    mol = Chem.MolFromSmiles(
        molecule.get_representation(MolecularGraph).canonical_smiles
    )

    fingerprints = list(
        AllChem.GetMorganGenerator(radius=radius)
        .GetSparseCountFingerprint(mol)
        .GetNonzeroElements()
        .keys()
    )
    print(fingerprints)
    return fingerprints


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
        name: str = "chembl",
    ):
        """
        Init FilterECFP with the path to the reference database and the radius.

        Use file "external_data/ecfp4_ChEMBL_ZINC.txt" for ZINC database.

        Args:
            path_db (str): file that contains the reference list of ECFC4
                fingerprints keys.
            radius (int, optional): radius if the ECFP fingerprint
                (2 for ECFP4, 1 for ECFP2). Defaults to 2.
            name (str, optional): name of the evaluation. Defaults to "chembl".
                use "chembl_zinc" for ZINC database.
        """
        super().__init__(f"UnknownECFP_{name}")

        self.radius = radius

        with open(path_db, encoding="utf8") as f:
            self.ecfp_list: frozenset[int] = frozenset(
                int(line.strip()) for line in f.readlines()
            )

        # create the Morgan fingerprints generator
        self.fingerprints_generator = AllChem.GetMorganGenerator(radius=2)

    @override
    def _evaluate(self, molecule: Molecule) -> int:

        mol = Chem.MolFromSmiles(
            molecule.get_representation(MolecularGraph).canonical_smiles
        )

        # get the Morgan fingerprints
        fingerprints = self.fingerprints_generator.GetSparseCountFingerprint(
            mol
        ).GetNonzeroElements()

        # count the number of ecfp that are not in the reference database
        nb_unknown_ecfp = 0
        for ecfp in fingerprints:
            if ecfp not in self.ecfp_list:
                nb_unknown_ecfp += 1

        return nb_unknown_ecfp


class FilterUnknownECFP(Evaluation):
    """Filter molecules with too many unknown ECFP."""

    def __init__(self, threshold: int = 0, name: str = "chembl"):
        super().__init__("FilterUnknownECFP")
        self.eval_name: str = f"UnknownECFP_{name}"
        self.threshold = threshold

    @override
    def _evaluate(self, molecule: Molecule) -> bool:
        try:
            value = molecule.value(self.eval_name)
        except KeyError as e:
            raise KeyError(
                "The molecule does not have a UnknownECFP value."
                "Make sure to evaluate it first with the UnknownECFP evaluation."
            ) from e

        if value <= self.threshold:
            return True

        raise EvaluationError(
            "The molecule contains more unknown ECFP than allowed "
            f"({value} > {self.threshold})."
        )
