"""
Functions to set default parameters for evomol.
"""

import os
from typing import Optional

from evomol import evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.representation import SMILES, MolecularGraph, Molecule


def setup_default_parameters(
    accepted_atoms: Optional[list[str]] = None, max_heavy_atoms: int = 38
) -> None:
    """Set the default parameters for evomol.

    Args:
        accepted_atoms (Optional[list[str]], optional): List of accepted atoms
            for the molecules. Defaults to None (["C", "O", "N", "F", "S"]).
        max_heavy_atoms (int, optional): Maximum number of heavy atoms in a
            molecule. Defaults to 38.
    """

    if accepted_atoms is None:
        accepted_atoms = ["C", "O", "N", "F", "S"]

    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = max_heavy_atoms

    Molecule.accepted_atoms = accepted_atoms


def setup_default_action_space(
    with_add_group: bool = False, with_remove_group: bool = False
) -> None:
    """Set the default action space for evomol."""

    mg.AddGroupMG.groups = [
        mg.Group("C1=CC=CS1", 5, [0]),
        mg.Group("C1=CC=CC=C1", 6, [0]),
        mg.Group("[N+](=O)[O-]", 3, [0]),
        mg.Group("N=[N+]=[N-]", 3, [0]),
        mg.Group("S(=O)(=O)O", 4, [0]),
    ]

    mg.ChangeBondMG.avoid_bond_breaking = False

    mg.RemoveGroupMG.remove_only_smallest = True
    mg.RemoveGroupMG.remove_only_single_bond = False
    mg.RemoveGroupMG.remove_charged_or_radical = True

    MolecularGraph.action_space = [
        mg.AddAtomMG,
        mg.ChangeBondMG,
        mg.CutAtomMG,
        mg.InsertCarbonMG,
        mg.MoveGroupMG,
        mg.RemoveAtomMG,
        mg.SubstituteAtomMG,
    ]
    if with_add_group:
        MolecularGraph.action_space.append(mg.AddGroupMG)
    if with_remove_group:
        MolecularGraph.action_space.append(mg.RemoveGroupMG)


def setup_filters(filter_name: str) -> list[evaluator.Evaluation]:
    """Return the evaluations for the generic cyclic features and ECFP

    Args:
        filter_name (str): chembl or chembl_zinc as filter

    Returns:
        list[evaluator.Evaluation]: list of evaluations
    """
    if filter_name == "chembl":
        return [
            evaluator.UnknownGCF(
                path_db=os.path.join("external_data", "gcf1.txt"),
                name="chembl",
            ),
            evaluator.FilterUnknownGCF(
                threshold=0,
                name="chembl",
            ),
            evaluator.UnknownECFP(
                path_db=os.path.join("external_data", "ecfp4_ChEMBL.txt"),
                radius=2,
                name="chembl",
            ),
            evaluator.FilterUnknownECFP(
                threshold=0,
                name="chembl",
            ),
        ]
    if filter_name == "chembl_zinc":
        return [
            evaluator.UnknownGCF(
                path_db=os.path.join("external_data", "gcf2.txt"),
                name="chembl_zinc",
            ),
            evaluator.FilterUnknownGCF(
                threshold=0,
                name="chembl_zinc",
            ),
            evaluator.UnknownECFP(
                path_db=os.path.join("external_data", "ecfp4_ChEMBL_ZINC.txt"),
                radius=2,
                name="chembl_zinc",
            ),
            evaluator.FilterUnknownECFP(
                threshold=0,
                name="chembl_zinc",
            ),
        ]
    raise ValueError("filter_name must be 'chembl' or 'chembl_zinc'")


def chembl_and_chembl_zinc_evaluations() -> list[evaluator.Evaluation]:
    """Return the evaluations for the generic cyclic features and ECFP

    Returns:
        list[evaluator.Evaluation]: list of evaluations
    """
    return [
        evaluator.UnknownGCF(
            path_db=os.path.join("external_data", "gcf1.txt"),
            name="chembl",
        ),
        evaluator.UnknownECFP(
            path_db=os.path.join("external_data", "ecfp4_ChEMBL.txt"),
            radius=2,
            name="chembl",
        ),
        evaluator.UnknownGCF(
            path_db=os.path.join("external_data", "gcf2.txt"),
            name="chembl_zinc",
        ),
        evaluator.UnknownECFP(
            path_db=os.path.join("external_data", "ecfp4_ChEMBL_ZINC.txt"),
            radius=2,
            name="chembl_zinc",
        ),
    ]
