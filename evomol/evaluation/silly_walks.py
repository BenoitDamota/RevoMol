"""
Silly Walks evaluation function.
"""

import json
from typing import Any

from rdkit import Chem
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation import MolecularGraph, Molecule


class SillyWalks(Evaluation):
    """
    Counting the proportion of bits in the ECFP4 fingerprint that never appear
    in the ChemBL.
    Based on the work of Patrick Walters (https://github.com/PatWalters/silly_walks)

    MIT License

    Copyright (c) 2020 Patrick Walters

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """

    def __init__(
        self,
        path_db: str = "external_data/evo_filters_data/complete_ChEMBL_ecfp4_dict.json",
        radius: int = 2,
    ):
        """
        Init SillyWalks with the path to the reference database and the radius.

        "external_data/evo_filters_data/complete_ChEMBL_ZINC_union_ecfp4_dict.
        json"
        if ZINC and CHEMBL

        Args:
            path_db (str): file that contains the reference dictionary of ECFC4
                           fingerprints keys
            radius (int, optional): radius if the ECFP fingerprint
                                    (2 for ECFP4, 1 for ECFP2).
                                    Defaults to 2.
        """
        super().__init__("SillyWalks")

        self.radius = radius

        # Reading the reference database
        with open(path_db, encoding="utf8") as f:
            self.count_dict: dict[str, list[str]] = json.load(f)

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:

        mol = Chem.MolFromSmiles(
            molecule.get_representation(MolecularGraph).canonical_smiles
        )

        if mol is None:
            return {
                "score": 1,
            }

        fingerprint = Chem.AllChem.GetMorganFingerprint(mol, self.radius)
        on_bits = fingerprint.GetNonzeroElements().keys()

        silly_bits = [
            x for x in [self.count_dict.get(str(ecfp)) for ecfp in on_bits] if not x
        ]
        score = len(silly_bits) / len(on_bits) if len(on_bits) > 0 else 0

        return {
            "score": score,
        }
