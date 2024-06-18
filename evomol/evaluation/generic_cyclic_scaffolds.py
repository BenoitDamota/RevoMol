"""
Evaluation of the number of unknown generic cyclic scaffolds in a molecule
compared to a reference dataset of generic cyclic scaffolds.
"""

from typing import Any

from rdkit import Chem
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation
from evomol.representation.molecular_graph import MolecularGraph
from evomol.representation.molecule import Molecule


def compute_generic_scaffold(smiles: str) -> str:
    """Computes the generic scaffold of a molecule.

    Args:
        smiles (str): SMILES representation of the molecule

    Returns:
        str: SMILES representation of the generic scaffold
    """
    # compute generic scaffold
    mol = Chem.MolFromSmiles(smiles)
    smiles_scaffolds: str = ""
    try:
        gscaf = Chem.Scaffolds.MurckoScaffold.MakeScaffoldGeneric(
            Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
        )
        smiles_scaffolds = Chem.MolToSmiles(gscaf)
    except Exception:  # pylint: disable=broad-exception-caught
        smiles_scaffolds = Chem.Scaffolds.MurckoScaffold.MurckoScaffoldSmiles(
            mol=Chem.MolFromSmiles(smiles), includeChirality=False
        )
        mol_scaffolds = MolecularGraph(Chem.MolFromSmiles(smiles_scaffolds))

        mol_scaffolds.update_representation()

        # Setting all bonds to single
        for atom1, atom2 in mol_scaffolds.bonds:
            mol_scaffolds.set_bond(atom1, atom2, 1)

        atom: Chem.Atom
        for atom in mol_scaffolds.mol.GetAtoms():
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
        smiles = Chem.MolToSmiles(gscaf)
    except Exception:  # pylint: disable=broad-exception-caught
        pass

    mol_graph = MolecularGraph(smiles)
    mol_graph.update_representation()

    # delete all bridge bonds
    for atom1, atom2 in mol_graph.bridge_bonds:
        mol_graph.set_bond(atom1, atom2, 0)

    mol_graph.update_representation()

    # computes smiles of rings features and remove unique atoms
    features = set(s for s in mol_graph.canonical_smiles.split(".") if len(s) > 4)

    return set(compute_generic_scaffold(feature) for feature in features)


class GenericCyclicScaffolds(Evaluation):
    """
    Returning the number of unknown generic cyclic scaffolds in given molecule
    compared to a reference dataset of generic cyclic scaffolds.

    Generic cyclic scaffolds are all rings extracted from the molecule and
    then converted so that
    * Each atom becomes a carbon atom
    * Each bond becomes a single bond
    """

    def __init__(self, path_db: str = "external_data/gcf2.txt"):
        """
        "external_data/gcf1.txt" if CHEMBL and not ZINC
        """
        super().__init__("GenericCyclicScaffolds")

        with open(path_db, encoding="utf8") as f:
            self.smiles_list: frozenset[str] = frozenset(
                line.strip() for line in f.readlines()
            )

    @override
    def _evaluate(self, molecule: Molecule) -> dict[str, Any]:

        mol_graph = molecule.get_representation(MolecularGraph)

        unique_features = set(extract_generic_scaffolds(mol_graph.canonical_smiles))

        unknown_count = sum(
            1 for feature in unique_features if feature not in self.smiles_list
        )

        return {
            "count": unknown_count,
        }
