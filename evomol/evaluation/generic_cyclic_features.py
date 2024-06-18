"""
This module contains the implementation of the GenericCyclicFeatures evaluation.
"""

from typing import Any

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation import MolecularGraph, Molecule


def compute_generic_scaffold_with_insat(smiles: str) -> str:
    """
    Compute a generic cyclic feature (with insaturations) of a cyclic subgraph.

    For that, it computes the Murcko scaffold on a subgraph (that includes a
    cycle), then all atoms with a coordination number of 4 or less are converted
    to carbon atoms. Since hypervalent carbon produce RDkit errors, hypervalent
    atoms are left unchanged.
    """
    smiles_scaf: str = MurckoScaffold.MurckoScaffoldSmiles(
        mol=Chem.MolFromSmiles(smiles), includeChirality=False
    )
    mol_scaf: MolecularGraph = MolecularGraph(smiles_scaf)

    for atom in mol_scaf.mol.GetAtoms():
        if atom.GetTotalValence() <= 4:
            # Changing atomic number
            atom.SetAtomicNum(6)
            # Setting formal charge to 0
            atom.SetFormalCharge(0)
        # else keep the atom type (eg. P with valence = 5)

    mol_scaf.update_representation()

    return mol_scaf.smiles


def compute_generic_cyclic_features_with_insat(smiles: str) -> list[str]:
    """
    Compute all generic cyclic features (with insaturations) of a molecules.

    For that, it breaks all the vertices (bonds) that do not belong to a cycle
    are deleted (called bridge in graph vocabulary).
    Then, compute_generic_scaffold_with_insat is called on all remaining cyclic
    subgraphs.

    """
    mol_graph = MolecularGraph(smiles)

    for atom1, atom2 in mol_graph.bridge_bonds:
        mol_graph.set_bond(atom1, atom2, 0)

    mol_graph.update_representation()

    smiles_subgraphs: str = Chem.MolToSmiles(mol_graph.mol).split(".")

    # computes smiles of rings features and remove unique atoms
    features: list[str] = [s for s in smiles_subgraphs if len(s) > 4]
    smiles_features: list[str] = [
        compute_generic_scaffold_with_insat(f) for f in features
    ]

    return smiles_features


class GenericCyclicFeatures(Evaluation):
    """https://github.com/BenoitDamota/gcf"""

    def __init__(self, path_db: str = "external_data/gcf2.txt"):
        """
        "external_data/gcf1.txt" if CHEMBL and not ZINC
        """
        super().__init__("GenericCyclicFeatures")

        with open(path_db, encoding="utf8") as f:
            self.smiles_list: frozenset[str] = frozenset(
                line.strip() for line in f.readlines()
            )

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:
        mol_graph = molecule.get_representation(MolecularGraph)

        unique_features = set(
            compute_generic_cyclic_features_with_insat(mol_graph.canonical_smiles)
        )

        unknown_count = sum(
            1 for feature in unique_features if feature not in self.smiles_list
        )

        return {
            "count": unknown_count,
        }
