"""
CLScore evaluation of a molecule using subgraph reference data.

Bühlmann, Sven, et Jean-Louis Reymond.
ChEMBL-Likeness Score and Database GDBChEMBL
Frontiers in Chemistry 8 (04/02/2020).
https://doi.org/10.3389/fchem.2020.00046.)

Based on https://github.com/reymond-group/GDBChEMBL
Options about, rooted, weighted, radius or cutoff have been removed and set to
default values.
"""

import os
import pickle

from rdkit import Chem
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation import MolecularGraph, Molecule


def extract_shingles(smiles: str) -> set[str]:
    """Extract shingles from a molecular graph.

    Args:
        smiles (str): canonical SMILES representation of the molecule.

    Returns:
        set[str]: The set of shingles.
    """
    shingles: set[str] = set()

    # Reloading molecule to make it aromatic
    mol: Chem.rdchem.RWMol = Chem.MolFromSmiles(smiles)

    radius = 3
    for atom in range(mol.GetNumAtoms()):
        for rad in range(1, radius + 1):
            bonds = Chem.FindAtomEnvironmentOfRadiusN(mol, rad, atom)

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
                    atom,
                    True,
                    False,
                    False,
                )
            )
    return shingles


class CLScore(Evaluation):
    """CLScore evaluation of a molecule using subgraph reference data."""

    def __init__(
        self,
        path: str = os.path.join(
            "external_data",
            "chembl_24_1_shingle_scores_log10_rooted_nchir_min_freq_100.pkl",
        ),
    ) -> None:
        super().__init__("CLscore")

        self.db_shingles: dict[str, float] = {}

        with open(path, "rb") as pyc:
            self.db_shingles = pickle.load(pyc)

    @override
    def _evaluate(self, molecule: Molecule) -> float:

        shingles = extract_shingles(
            molecule.get_representation(MolecularGraph).aromatic_canonical_smiles
        )

        if not shingles:
            return 0.0

        # using log10 of shingle frequency
        # if key not present, add 0 per default
        scores: float = sum(self.db_shingles.get(shingle, 0.0) for shingle in shingles)

        cl_score: float = scores / len(shingles)

        return cl_score
