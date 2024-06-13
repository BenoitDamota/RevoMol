import os
import pickle
from typing import Any

from rdkit import Chem
from typing_extensions import override

from evomol.evaluation import Evaluation
from evomol.representation.molecule import Molecule
from evomol.representation.molecular_graph import MolecularGraph


class CLScore(Evaluation):
    """Code from https://github.com/reymond-group/GDBChEMBL"""

    def __init__(
        self,
        path=os.path.join(
            "external_data",
            "chembl_24_1_shingle_scores_log10_rooted_nchir_min_freq_100.pkl",
        ),
        radius=3,
    ):
        super().__init__("CL_score")

        self.radius = radius
        self.rooted = True
        self.weighted = True
        self.cut_off = 0.0
        self.db_shingles: dict = {}

        # Loading ChEMBL shingles database
        # file: str = ""

        # if self.rooted:
        #     file = "chembl_24_1_shingle_scores_log10_rooted_nchir_min_freq_100.pkl"
        # else:
        #     file = "chembl_24_1_shingle_scores_log10_nrooted_nchir.pkl"
        # path = os.path.join(os.environ["SHINGLE_LIBS"], file)
        with open(path, "rb") as pyc:
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
            cl_score = 0

        return {
            "CL_score": cl_score,
        }
