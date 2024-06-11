from abc import ABC, abstractmethod
import os
import sys
import pickle
import json
from typing import Any

from rdkit import Chem, AllChem
import networkx as nx
import guacamol
from typing_extensions import override
import pandas as pd

from evomol.representation.molecule import Molecule
from evomol.representation.molecular_graph import MolecularGraph

sys.path.append(os.path.join(Chem.RDConfig.RDContribDir, "SA_Score"))
import sascorer  # pylint: disable=import-error, wrong-import-position, wrong-import-order


class Evaluation(ABC):

    def __init__(self, name: str):
        self._name = name

    @property
    def name(self):
        return self._name

    def evaluate(self, molecule: Molecule) -> dict[str, Any]:
        value = molecule.value(self.name)
        if value is not None:
            return value
        value = self._evaluate(molecule)
        molecule.set_value(self.name, value)
        return value

    @abstractmethod
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        pass


class LogP(Evaluation):

    def __init__(self):
        super().__init__("logP")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        # normalization constants
        # statistics from 250k_rndm_zinc_drugs_clean.smi
        logP_mean = 2.4570953396190123
        logP_std = 1.434324401111988

        mol_graph = molecule.get_representation(MolecularGraph)

        log_p = Chem.Descriptors.MolLogP(mol_graph)

        normalized_log_p = (log_p - logP_mean) / logP_std

        return {
            "logP": log_p,
            "normalized_logP": normalized_log_p,
        }


class SAScore(Evaluation):

    def __init__(self):
        super().__init__("SA_score")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        mol_graph = molecule.get_representation(MolecularGraph)

        sa_score = -sascorer.calculateScore(mol_graph)

        # normalization constants
        # statistics from 250k_rndm_zinc_drugs_clean.smi
        sa_mean = -3.0525811293166134
        sa_std = 0.8335207024513095

        normalized_sa_score = (sa_score - sa_mean) / sa_std

        # Ertl, Peter, et Ansgar Schuffenhauer.
        # Estimation of synthetic accessibility score of drug-like
        # molecules based on molecular complexity and fragment contributions
        # Journal of Cheminformatics 1, no 1 (10 juin 2009): 8.
        # https://doi.org/10.1186/1758-2946-1-8.
        other_sa_score = 1 - (sa_score - 1) / 9

        return {
            "SA_score": sa_score,
            "normalized_SA_score": normalized_sa_score,
            "other_SA_score": other_sa_score,
        }


class CycleScore(Evaluation):

    def __init__(self):
        super().__init__("Cycle_score")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        mol_graph = molecule.get_representation(MolecularGraph)

        # normalization constants
        # statistics from 250k_rndm_zinc_drugs_clean.smi
        cycle_mean = -0.0485696876403053
        cycle_std = 0.2860212110245455

        # cycle score
        cycles = nx.cycle_basis(nx.Graph(mol_graph.adjacency_matrix))
        # or
        cycles = mol_graph.mol.GetRingInfo().AtomRings()

        cycle_length = max(len(cycle) for cycle in cycles) if cycles else 0

        cycle_length = max(cycle_length - 6, 0)

        cycle_score = -cycle_length

        normalized_cycle = (cycle_score - cycle_mean) / cycle_std

        return {
            "cycle_score": cycle_score,
            "normalized_cycle_score": normalized_cycle,
        }


class CLScore(Evaluation):

    def __init__(self):
        super().__init__("CL_score")

        self.scores = None
        self.radius = 3
        self.rooted = True
        self.weighted = True
        self.cut_off = 0.0
        self.db_shingles: dict = {}

        # Loading ChEMBL shingles database
        file: str = ""

        if self.rooted:
            file = "chembl_24_1_shingle_scores_log10_rooted_nchir_min_freq_100.pkl"
        else:
            file = "chembl_24_1_shingle_scores_log10_nrooted_nchir.pkl"

        with open(os.path.join(os.environ["SHINGLE_LIBS"], file), "rb") as pyc:
            self.db_shingles = pickle.load(pyc)

    def extract_shingles(self, mol_graph: MolecularGraph):
        shingles = set()

        # Reloading molecule to make it aromatic
        mol: Chem.rdchem.RWMol = Chem.MolFromSmiles(mol_graph.canonical_smiles)

        for atom in range(mol.GetNumAtoms()):
            for radius in range(1, self.radius + 1):
                bonds = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom)

                if not bonds:
                    break

                atoms = set()
                for bond_idx in bonds:
                    bond = mol.GetBondWithIdx(bond_idx)
                    atoms.add(bond.GetBeginAtomIdx())
                    atoms.add(bond.GetEndAtomIdx())

                shingles.add(
                    Chem.rdmolfiles.MolFragmentToSmiles(
                        mol,
                        list(atoms),
                        bonds,
                        0,
                        0,
                        False,
                        False,
                        atom if self.rooted else -1,
                        True,
                        False,
                        False,
                    )
                )
        return shingles

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:

        shingles = self.extract_shingles(molecule.get_representation(MolecularGraph))

        if not shingles:
            return {
                "CL_score": 0,
            }

        scores = 0
        if self.weighted:
            # using log10 of shingle frequency
            # if key not present, add 0 per default
            scores = sum(self.db_shingles.get(shingle, 0) for shingle in shingles)
        else:
            # working binary (i.e. if present -> count++ )
            scores = sum(1 for shingle in shingles if shingle in self.db_shingles)

        cl_score = scores / len(shingles)

        if not (self.cut_off == 0.0 or self.cut_off <= cl_score):
            cl_score = None

        return {
            "CL_score": cl_score,
        }


class IsomerGuacamol(Evaluation):
    """
    Isomer score based on the implementation of GuacaMol
    Nathan Brown et al.
    GuacaMol: Benchmarking Models for de Novo Molecular Design
    Journal of Chemical Information and Modeling 59, no. 3
    (March 25, 2019): 1096–1108
    https://doi.org/10.1021/acs.jcim.8b00839
    """

    def __init__(self, formula: str):
        super().__init__("isomer_guacamol")
        self.scorer = guacamol.common_scoring_functions.IsomerScoringFunction(formula)

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        mol_graph = molecule.get_representation(MolecularGraph)

        return {
            "isomer_score": self.scorer.score(mol_graph.canonical_smiles),
        }


class Rediscovery(Evaluation):
    """
    Rediscovery score based on the implementation of GuacaMol
    Nathan Brown et al.
    GuacaMol: Benchmarking Models for de Novo Molecular Design
    Journal of Chemical Information and Modeling 59, no. 3 (March 25, 2019): 1096–1108
    https://doi.org/10.1021/acs.jcim.8b00839
    """

    def __init__(self, target_smiles: str):
        super().__init__("rediscovery_" + target_smiles)
        self.scorer = guacamol.common_scoring_functions.TanimotoScoringFunction(
            target_smiles,
            fp_type="ECFP4",
        )

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        mol_graph = molecule.get_representation(MolecularGraph)

        return {
            "rediscovery_score": self.scorer.score(mol_graph.canonical_smiles),
        }


class QED(Evaluation):
    """
    Evaluation of population with QED score using RDKit implementation.
    Bickerton, G. Richard, Gaia V. Paolini, Jérémy Besnard, Sorel Muresan, et Andrew L. Hopkins.
    Quantifying the Chemical Beauty of Drugs
    Nature Chemistry 4, nᵒ 2 (février 2012): 90-98
    https://doi.org/10.1038/nchem.1243
    """

    def __init__(self):
        super().__init__("QED")

    @override
    def _evaluate(self, molecule: Molecule) -> None:

        mol_graph: MolecularGraph = molecule.get_representation(MolecularGraph)

        score = Chem.QED.qed(Chem.MolFromSmiles(mol_graph.canonical_smiles))

        return {
            "QED_score": score,
        }


def compute_generic_scaffold(smiles: str) -> str:
    # compute generic scaffold
    mol = Chem.MolFromSmiles(smiles)
    try:
        gscaf = Chem.Scaffolds.MurckoScaffold.MakeScaffoldGeneric(
            Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
        )
        smiles_scaffolds = Chem.MolToSmiles(gscaf)
    except Exception:
        smiles_scaffolds = Chem.Scaffolds.MurckoScaffold.MurckoScaffoldSmiles(
            mol=Chem.MolFromSmiles(smiles), includeChirality=False
        )
        mol_scaffolds = MolecularGraph(Chem.MolFromSmiles(smiles_scaffolds))

        mol_scaffolds.update_representation()

        # Setting all bonds to single
        for atom1, atom2 in mol_scaffolds.bonds:
            mol_scaffolds.set_bond(atom1, atom2, 1)

        for atom in mol_scaffolds.mol_graph.GetAtoms():
            if atom.GetTotalValence() <= 4:
                # Changing atomic number
                atom.SetAtomicNum(6)
                # Setting formal charge to 0
                atom.SetFormalCharge(0)

        mol_scaffolds.update_representation()
        smiles_scaffolds = mol_scaffolds.canonical_smiles
    return smiles_scaffolds


def extract_generic_scaffolds(smiles: str) -> set[str]:
    """
    Returning the list of generic cyclic scaffolds for given SMILES.

    Generic cyclic scaffolds are all rings extracted from the molecule and then
    converted so that
    * Each atom becomes a carbon atom
    * Each bond becomes a single bond

    If an atom shares 4 cyclic neighbors or more,
    its type is kept unchanged since it cannot be converted to a carbon atom.
    """
    mol = Chem.MolFromSmiles(smiles)
    try:
        gscaf = Chem.Scaffolds.MurckoScaffold.MakeScaffoldGeneric(
            Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
        )
        mol = Chem.MolToSmiles(gscaf)
    except Exception:
        pass

    mol_graph = MolecularGraph()
    mol_graph.mol = mol
    mol_graph.update_representation()

    # delete all bridge bonds
    for atom1, atom2 in mol_graph.bridge_bonds:
        mol_graph.set_bond(atom1, atom2, 0)

    mol_graph.update_representation()

    # computes smiles of rings features and remove unique atoms
    features = set(s for s in mol_graph.canonical_smiles.split(".") if len(s) > 4)

    return set(compute_generic_scaffold(feature) for feature in features)


class UnkownGenericCyclicScaffolds(Evaluation):
    """
    Returning the number of unknown generic cyclic scaffolds in given molecule
    compared to a reference dataset of generic cyclic scaffolds.

    Generic cyclic scaffolds are all rings extracted from the molecule and
    then converted so that
    * Each atom becomes a carbon atom
    * Each bond becomes a single bond
    """

    def __init__(self, path_db: str):
        super().__init__("UnkownGenericCyclicScaffolds")

        with open(path_db, "r", encoding="utf8") as f:
            self.ref_dic: dict = json.load(f)

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:

        mol_graph = molecule.get_representation(MolecularGraph)

        unique_features = set(extract_generic_scaffolds(mol_graph.canonical_smiles))

        unknown_count = sum(
            1 for feature in unique_features if feature not in self.ref_dic
        )

        return {
            "count": unknown_count,
        }


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

    def __init__(self, path_db: str, radius: int = 2):
        """
        Init SillyWalks with the path to the reference database and the radius.

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
        with open(path_db, "r", encoding="utf8") as f:
            self.count_dic: dict = json.load(f)

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:

        mol = Chem.MolFromSmiles(
            molecule.get_representation(MolecularGraph).canonical_smiles
        )

        if mol is None:
            return {
                "score": 1,
            }

        fingerprint = AllChem.GetMorganFingerprint(mol, self.radius)
        on_bits = fingerprint.GetNonzeroElements().keys()

        silly_bits = [x for y in on_bits for x in self.count_dict.get(str(y)) if not x]
        score = len(silly_bits) / len(on_bits) if len(on_bits) > 0 else 0

        return {
            "score": score,
        }


class NPerturbations(Evaluation):
    """
    Strategy that computes the count of perturbations that were already applied
    to the solutions.
    If the parameters of EvoMol are set such as a mutation is a single
    perturbation on the molecular graph, then this value is equivalent to the
    number of previous mutations.
    """

    def __init__(self):
        super().__init__("NPerturbations")

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        return {
            "count": molecule.nb_perturbations,
        }


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

    def __init__(self):
        super().__init__("RDFilters")

        rules_file_name: str = os.environ["FILTER_RULES_DATA"] + "/rules.json"

        alert_file_name: str = os.environ["FILTER_RULES_DATA"] + "/alert_collection.csv"

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
