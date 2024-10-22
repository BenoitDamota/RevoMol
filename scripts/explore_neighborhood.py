"""
This script is used to print the neighborhood of a molecule at a depth of 1.
It also shows the evaluation of the molecule with the GCF and ECFP evaluators.
"""

import os
import sys

import typer

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol import evaluation as evaluator
from evomol.action import pprint_action_space
from evomol.representation import MolecularGraph, Molecule


def explore_neighborhood(smiles: str) -> None:
    """Explore the neighborhood of a molecule at a depth of 1.

    Args:
        smiles (str): SMILES representation of the molecule.
    """
    # init default parameters and action space
    dp.setup_default_parameters()
    dp.setup_default_action_space(with_add_group=False, with_remove_group=False)

    # Load the GCF and ECFP evaluators
    gcf_1 = evaluator.UnknownGCF(
        path_db=os.path.join("external_data", "gcf1.txt"), name="chembl"
    )
    ecfp_1 = evaluator.UnknownECFP(
        path_db=os.path.join("external_data", "ecfp4_ChEMBL.txt"),
        radius=2,
        name="chembl",
    )
    gcf_2 = evaluator.UnknownGCF(
        path_db=os.path.join("external_data", "gcf2.txt"), name="chembl_zinc"
    )
    ecfp_2 = evaluator.UnknownECFP(
        path_db=os.path.join("external_data", "ecfp4_ChEMBL_ZINC.txt"),
        radius=2,
        name="chembl_zinc",
    )

    mol = Molecule(smiles)
    can_smiles = mol.get_representation(MolecularGraph).canonical_smiles
    print(f"{mol=} {can_smiles=}")

    # compute the possible actions and print them
    possible_actions = mol.compute_possible_actions()
    pprint_action_space(possible_actions)

    legal_chembl = []
    legal_chembl_zinc = []
    neighborhood_list: list[str] = []
    # apply each action and print the new molecule
    for _, action_space in possible_actions.items():
        for _, actions_list in action_space.items():
            for action in actions_list:
                print(action, "->", end=" ")
                new_mol = action.apply()
                print(new_mol, end=f"{'':4}")
                neighborhood_list.append(
                    new_mol.get_representation(MolecularGraph).canonical_smiles
                )

                # evaluate the new molecule with the GCF and ECFP evaluators
                nb_gcf1 = gcf_1.evaluate(new_mol)
                nb_gcf2 = gcf_2.evaluate(new_mol)
                nb_ecfp1 = ecfp_1.evaluate(new_mol)
                nb_ecfp2 = ecfp_2.evaluate(new_mol)
                print(f"{nb_gcf1=} {nb_gcf2=} {nb_ecfp1=} {nb_ecfp2=}")
                if nb_gcf1 == 0 and nb_ecfp1 == 0:
                    legal_chembl.append(new_mol)
                if nb_gcf2 == 0 and nb_ecfp2 == 0:
                    legal_chembl_zinc.append(new_mol)

    # print the number of legal molecules
    print(f"Legal molecules in chembl: {len(legal_chembl)}")
    print(f"Legal molecules in chembl set: {len(set(legal_chembl))}")
    print(f"Legal molecules in chembl_zinc: {len(legal_chembl_zinc)}")
    print(f"Legal molecules in chembl_zinc set: {len(set(legal_chembl_zinc))}")
    print(f"Neighborhood size: {len(neighborhood_list)}")
    print(f"Neighborhood set: {len(set(neighborhood_list))}")


if __name__ == "__main__":
    # smi = "C"
    # smi = "C1=C2NNS13C=C23"
    # explore_neighborhood(smi)

    typer.run(explore_neighborhood)
