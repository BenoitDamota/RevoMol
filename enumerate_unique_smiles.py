import time
import sys

import typer

from evomol import evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.representation import SMILES, MolecularGraph, Molecule
from evomol.action import Action


def is_valid_molecule(
    molecule: Molecule, evaluations: list[evaluator.Evaluation]
) -> bool:
    """Check if a molecule is valid

    Args:
        molecule (Molecule): Molecule to check
        evaluations (list[evaluator.Evaluation]): filter to apply

    Returns:
        bool: True if the molecule is valid
    """
    for eval_ in evaluations:
        try:
            eval_.evaluate(molecule)
        except evaluator.EvaluationError:
            return False
    return True


def list_neighbors(molecule: Molecule) -> list[Molecule]:
    new_mols = []
    molecule.list_all_possible_actions()
    while molecule.nb_possible_actions() != 0:
        action_space = list(molecule.possible_actions.keys())[0]
        while (
            molecule.possible_actions is not None
            and molecule.possible_actions.get(action_space) is not None
        ):
            action = list(molecule.possible_actions[action_space].keys())[0]
            new_mols.append(molecule.possible_actions[action_space][action][0].apply())
    return new_mols


def find_unique_neighbors(
    molecule: Molecule,
    depth: int,
    found_valid: dict[Molecule, int],
    found_invalid: dict[Molecule, int],
):
    to_explore: set[Molecule] = {molecule}
    explored: set[Molecule] = set()

    current_depth = 1

    while current_depth <= depth:
        total_new_mols = []
        while to_explore:
            current_mol = to_explore.pop()

            total_new_mols.extend(list_neighbors(current_mol))

        for new_mol in total_new_mols:
            if new_mol in found_valid:
                found_valid[new_mol] += 1
                continue
            if new_mol in found_invalid:
                found_invalid[new_mol] += 1
            if new_mol in explored:
                continue

            explored.add(new_mol)
            if current_depth < depth:
                to_explore.add(new_mol)

        current_depth += 1

    return explored


"""
A -> B -> C -> D

Si on veut une profondeur de 3 de A:
    on va lister tous les voisins de niveau 1 (B,C)
    puis tous les voisins de niveau 2 (D,E,C)
    puis tous les voisins de niveau 3 (F,G,H,C,D,B)
    
    puis on va filtrer l'ensemble

"""


def explore_smiles_unique(
    starting_smiles: str, depth: int, evaluations: list[evaluator.Evaluation]
) -> dict[str, int]:
    """Explore the neighbors of a SMILES

    Args:
        starting_smiles (str): SMILES to explore
        depth (int): Maximum depth to explore
        evaluations (list[evaluator.Evaluation]): filter to apply

    Returns:
        dict[str, int]: SMILES and number of times it was found
    """

    molecule = Molecule(starting_smiles)

    # set of SMILES to explore
    to_explore: set[Molecule] = set([molecule])

    # dictionary of SMILES and number of times it was found
    found_valid: dict[Molecule, int] = {molecule: 1}
    found_invalid: dict[Molecule, int] = {}

    while to_explore:
        # get the next SMILES to explore
        current_molecule: Molecule = to_explore.pop()

        # find the neighbors of the current SMILES
        new_mols = find_unique_neighbors(
            current_molecule, depth, found_valid, found_invalid
        )

        # check if the neighbors are valid and add them to the set to_try
        for mol in new_mols:
            if mol in found_valid:
                continue
            if mol in found_invalid:
                continue

            if is_valid_molecule(mol, evaluations):
                to_explore.add(mol)
                found_valid[mol] = 1
            else:
                found_invalid[mol] = 1
    return found_valid


def enumerate_unique_from_smiles(
    start_smiles: str, eval_name: str, nb_heavy_atoms: int, depth: int
):
    """Enumerate the neighbors of a SMILES

    Args:
        start_smiles (str): SMILES to start from
        eval_name (str): chembl or chembl_zinc as filter
        nb_heavy_atoms (int): maximum number of heavy atoms
        depth (int): maximum depth to explore
    """
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.max_heavy_atoms = nb_heavy_atoms

    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    MolecularGraph.action_space = [
        mg.AddAtomMG,
        # mg.AddGroupMG,
        mg.ChangeBondMG,
        mg.CutAtomMG,
        mg.InsertCarbonMG,
        mg.MoveGroupMG,
        mg.RemoveAtomMG,
        # mg.RemoveGroupMG,
        mg.SubstituteAtomMG,
    ]

    evaluations = []
    if eval_name == "chembl":
        evaluations = [
            evaluator.UnknownGCF(path_db="external_data/gcf1.txt"),
            evaluator.FilterUnknownGCF(threshold=0),
            evaluator.UnknownECFP(path_db="external_data/ecfp4_ChEMBL.txt", radius=2),
            evaluator.FilterUnknownECFP(threshold=0),
        ]
    elif eval_name == "chembl_zinc":
        evaluations = [
            evaluator.UnknownGCF(path_db="external_data/gcf2.txt"),
            evaluator.FilterUnknownGCF(threshold=0),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL_ZINC.txt", radius=2
            ),
            evaluator.FilterUnknownECFP(threshold=0),
        ]

    start = time.time()

    found = explore_smiles_unique(start_smiles, depth, evaluations)

    duration = time.time() - start

    print("molecule,heavy_atom_limit,depth,eval,nb_molecules,time")
    print(
        ",".join(
            map(
                str,
                [
                    start_smiles,
                    nb_heavy_atoms,
                    depth,
                    eval_name,
                    len(found),
                    f"{duration:.2f}",
                ],
            )
        )
    )

    file = f"output/enumerate_unique_{start_smiles}_{eval_name}_{nb_heavy_atoms}.txt"

    with open(file, "w") as f:
        for mol, nb in found.items():
            if mol == "":
                f.write(f'"" {nb}\n')
            else:
                f.write(f"{mol} {nb}\n")


if __name__ == "__main__":
    # parallel -j 20 python enumerate_unique_smiles.py C {1} {2} 2 ::: chembl chembl_zinc ::: {1..10}
    print("this code doesn't work, there is some debugging to do")
    sys.exit(0)
    typer.run(enumerate_unique_from_smiles)
