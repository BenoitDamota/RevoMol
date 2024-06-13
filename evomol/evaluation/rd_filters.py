import os
from typing import Any
import json

from typing_extensions import override
from rdkit import Chem
import pandas as pd


from evomol.evaluation import Evaluation
from evomol.representation.molecule import Molecule
from evomol.representation.molecular_graph import MolecularGraph


class RDFilters(Evaluation):
    """
    Adapted from https://github.com/PatWalters/rd_filters

    MIT License

    Copyright (c) 2018 Patrick Walters

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    """

    def __init__(self, path="external_data"):
        super().__init__("RDFilters")

        rules_file_name: str = os.path.join(path, "rules.json")

        alert_file_name: str = os.path.join(path, "alert_collection.csv")

        # make sure there wasn't a blank line introduced
        rule_df = pd.read_csv(alert_file_name).dropna()

        with open(rules_file_name, encoding="utf8") as json_file:
            self.rule_dict = json.load(json_file)

        rules_list = [
            x.replace("Rule_", "")
            for x in self.rule_dict.keys()
            if x.startswith("Rule") and self.rule_dict[x]
        ]

        rule_df = rule_df[rule_df.rule_set_name.isin(rules_list)]

        tmp_rule_list = rule_df[
            ["rule_id", "smarts", "max", "description"]
        ].values.tolist()

        self.rule_list = []
        for _, smarts, max_val, desc in tmp_rule_list:
            smarts_mol = Chem.MolFromSmarts(smarts)
            if smarts_mol:
                self.rule_list.append([smarts_mol, max_val, desc])

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:

        mol: Chem.rdchem.RWMol = Chem.MolFromSmiles(
            molecule.get_representation(MolecularGraph).canonical_smiles
        )

        if mol is None:
            return {
                "filter": 0,
            }

        desc_list = [
            Chem.Descriptors.MolWt(mol),
            Chem.Descriptors.MolLogP(mol),
            Chem.Descriptors.NumHDonors(mol),
            Chem.Descriptors.NumHAcceptors(mol),
            Chem.Descriptors.TPSA(mol),
            Chem.rdMolDescriptors.CalcNumRotatableBonds(mol),
        ]
        df = pd.DataFrame(
            [desc_list], columns=["MW", "LogP", "HBD", "HBA", "TPSA", "Rot"]
        )
        df_ok = df[
            df.MW.between(*(self.rule_dict["MW"]))
            & df.LogP.between(*(self.rule_dict["LogP"]))
            & df.HBD.between(*(self.rule_dict["HBD"]))
            & df.HBA.between(*(self.rule_dict["HBA"]))
            & df.TPSA.between(*(self.rule_dict["TPSA"]))
        ]

        if len(df_ok) == 0:
            return {"filter": 0}

        for patt, max_val, _ in self.rule_list:
            if len(mol.GetSubstructMatches(patt)) > max_val:
                return {"filter": 0}

        return {"filter": 1}
