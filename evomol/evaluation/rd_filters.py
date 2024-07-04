"""
This module contains the RDFilters class, which is used to evaluate molecules
based on the RDKit filters.


Adapted from https://github.com/PatWalters/rd_filters
"""

import json
import os

import pandas as pd
from rdkit import Chem
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation, EvaluationError
from evomol.representation import MolecularGraph, Molecule


class RDFilters(Evaluation):
    """RDFilters evaluation class."""

    def __init__(self, path: str = "external_data") -> None:
        super().__init__("RDFilters")

        rules_file_name: str = os.path.join(path, "rules.json")

        alert_file_name: str = os.path.join(path, "alert_collection.csv")

        # make sure there wasn't a blank line introduced
        rule_df = pd.read_csv(alert_file_name).dropna()

        with open(rules_file_name, encoding="utf8") as json_file:
            self.rule_dict = json.load(json_file)

        rules_list = [
            x.replace("Rule_", "")
            for x in self.rule_dict
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
    def _evaluate(self, molecule: Molecule) -> bool:

        mol: Chem.rdchem.RWMol = Chem.MolFromSmiles(
            molecule.get_representation(MolecularGraph).aromatic_canonical_smiles
        )

        if mol is None:
            raise EvaluationError("RDFilters - The molecule is None.")

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
            raise EvaluationError("RDFilters - The molecule is out of range.")

        for patt, max_val, _ in self.rule_list:
            if len(mol.GetSubstructMatches(patt)) > max_val:
                raise EvaluationError(
                    f"RDFilters - The molecule has {patt} substructure."
                )

        return True
