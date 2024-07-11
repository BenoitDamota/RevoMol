"""
Evaluation of the number of unknown generic cyclic features in a molecule
compared to a reference dataset of generic cyclic features.

Adapted from https://github.com/BenoitDamota/gcf
"""

import os
from typing import Optional

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from typing_extensions import override

from evomol.evaluation.evaluation import Evaluation, EvaluationError
from evomol.representation import MolecularGraph, Molecule


def try_convert_to_carbon(smiles: str) -> Optional[str]:
    """Try to convert the atoms in the molecule to carbon atoms.
    Should not work if the molecule contains P or S atoms with valence > 4.

    Args:
        smiles (str): SMILES representation of the molecule

    Returns:
        Optional[str]: SMILES representation of the molecule with all atoms
            converted to carbon atoms or None if it fails
    """
    try:
        gscaf = MurckoScaffold.MakeScaffoldGeneric(
            MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(smiles))
        )

        smiles_scaffolds: str = Chem.MolToSmiles(gscaf)
        return smiles_scaffolds
    except Exception:  # pylint: disable=broad-exception-caught
        return None


def convert_to_carbon(smiles: str) -> str:
    """Convert the atoms in the molecule to carbon atoms, except for atoms with
    a valence > 4.

    Compute a generic cyclic feature (with insaturations) of a cyclic subgraph.

    For that, it computes the Murcko scaffold on a subgraph (that includes a
    cycle).
    Then all atoms with a coordination number of 4 or less are converted
    to carbon atoms.
    Since hypervalent carbon produce RDkit errors, hypervalent atoms are left
    unchanged.

    Args:
        smiles (str): SMILES representation of the molecule

    Returns:
        str: SMILES representation of the molecule with all atoms converted to
            carbon atoms
    """
    smiles_scaffolds = try_convert_to_carbon(smiles)
    if smiles_scaffolds is not None:
        return smiles_scaffolds

    smiles_scaffolds = MurckoScaffold.MurckoScaffoldSmiles(
        mol=Chem.MolFromSmiles(smiles), includeChirality=False
    )
    mol_scaffolds = MolecularGraph(smiles_scaffolds)

    atom: Chem.Atom
    for atom in mol_scaffolds.mol.GetAtoms():
        if atom.GetTotalValence() <= 4:
            # change atomic number
            atom.SetAtomicNum(6)
            # set formal charge to 0
            atom.SetFormalCharge(0)
        # else keep the atom type (eg. P with valence = 5)

    mol_scaffolds.update_representation()
    return mol_scaffolds.canonical_smiles


def list_gcf(smiles: str) -> set[str]:
    """Compute all generic cyclic features (with insaturations) of a molecule.

    All bridges (bonds that do not belong to a cycle) are removed.
    Then, compute_generic_scaffold_with_insat is called on all remaining cyclic
    subgraphs.

    Args:
        smiles (str): SMILES representation of the molecule

    Returns:
        list[str]: List of generic cyclic features (with insaturations)
    """
    mol_graph = MolecularGraph(smiles)

    # delete all bridge bonds
    for atom1, atom2 in mol_graph.bridge_bonds:
        mol_graph.set_bond(atom1, atom2, 0)

    mol_graph.update_representation()

    subgraphs: set[str] = set(mol_graph.canonical_smiles.split("."))

    # computes smiles of rings features and remove unique atoms
    # the length condition is to avoid computing the scaffold of subgraphs
    # that are not cyclic (minimum length of the SMILES of a cycle is 5,
    # as it is a cycle with 3 atoms and 2 numbers representing the ring closure
    # , eg. C1CC1)
    features: set[str] = set(
        convert_to_carbon(smiles) for smiles in subgraphs if len(smiles) > 4
    )

    return features


class UnknownGCF(Evaluation):
    """
    GenericCyclicFeatures"""

    def __init__(
        self,
        path_db: str = os.path.join("external_data", "gcf1.txt"),
        name: str = "chembl",
    ):
        """Init UnknownGCF with the path to the reference database.

        Args:
            path_db (str, optional): file that contains the reference list of
                generic cyclic features (GCF) keys.
                Defaults to os.path.join("external_data", "gcf1.txt").
                Use "external_data/gcf2.txt" for chembl_zinc database.
            name (str, optional): name of the database used.
                Defaults to "chembl" for gcf1.txt.
                Use "chembl_zinc" for gcf2.txt.
        """

        super().__init__(f"GenericCyclicFeatures_{name}")

        with open(path_db, encoding="utf8") as f:
            self.smiles_list: frozenset[str] = frozenset(
                line.strip() for line in f.readlines()
            )

    @override
    def _evaluate(self, molecule: Molecule) -> int:
        mol_graph = molecule.get_representation(MolecularGraph)

        unique_features = list_gcf(mol_graph.canonical_smiles)

        unknown_count = 0
        for feature in unique_features:
            if feature not in self.smiles_list:
                unknown_count += 1

        return unknown_count


class FilterUnknownGCF(Evaluation):
    """Filter molecules with too many unknown GenericCyclicFeatures."""

    def __init__(self, threshold: int = 0):
        super().__init__("FilterUnknownGCF")
        self.threshold = threshold

    @override
    def _evaluate(self, molecule: Molecule) -> bool:
        value: int = molecule.value("GenericCyclicFeatures")
        if value is None:
            raise RuntimeError(
                "The molecule does not have a GenericCyclicFeatures value."
            )
        if value > self.threshold:
            raise EvaluationError(
                "The molecule contains more unknown GCF than allowed "
                f"({self.threshold})."
            )
        return True
