import time

import typer

from evomol import evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.representation import SMILES, MolecularGraph, Molecule


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


def find_neighbors(molecule: Molecule, max_depth: int) -> set[Molecule]:
    """Iteratively find the neighbors of a molecule up to a certain depth

    Args:
        molecule (Molecule): Molecule to explore
        max_depth (int): Maximum depth to explore

    Returns:
        set[Molecule]: set of molecules found
    """
    neighbors: set = set()
    queue = {molecule}
    depth = 1

    while depth <= max_depth:
        next_queue = set()
        for current_mol in queue:
            for new_mol in list_neighbors(current_mol):
                if new_mol not in neighbors:
                    next_queue.add(new_mol)
                    neighbors.add(new_mol)
        queue = next_queue
        depth += 1

    return neighbors


def explore_smiles(
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

    # set of SMILES to explore
    to_explore: set[str] = {starting_smiles}

    # dictionaries of SMILES and number of times it was found
    found_valid: dict[str, int] = {starting_smiles: 1}
    found_invalid: dict[str, int] = {}

    while to_explore:
        # get the next SMILES to explore
        current_smiles: str = to_explore.pop()

        # find the neighbors of the current SMILES
        new_mols = find_neighbors(Molecule(current_smiles), depth)

        # check if the neighbors are valid and add them to the set to_try
        for new_mol in new_mols:
            new_smi = str(new_mol)
            if new_smi in found_valid:
                found_valid[new_smi] += 1
            elif new_smi in found_invalid:
                found_invalid[new_smi] += 1
            elif is_valid_molecule(new_mol, evaluations):
                to_explore.add(new_smi)
                found_valid[new_smi] = 1
            else:
                found_invalid[new_smi] = 1

    return found_valid


def enumerate_from_smiles(
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

    print("molecule,heavy_atom_limit,depth,eval,nb_molecules,time")

    start_time = time.time()

    can_smi = Molecule(start_smiles).get_representation(MolecularGraph).canonical_smiles

    found = explore_smiles(can_smi, depth, evaluations)

    duration = time.time() - start_time

    print(f"{can_smi},{nb_heavy_atoms},{depth},{eval_name},{len(found)},{duration:.2f}")

    file_path = f"output/enumerate_{can_smi}_{eval_name}_{nb_heavy_atoms}.txt"
    with open(file_path, "w") as file:
        for mol, nb in found.items():
            file.write(f"{mol} {nb}\n")


if __name__ == "__main__":
    # parallel -j 20 python enumerate_smiles.py C {1} {2} 2 ::: chembl chembl_zinc ::: {1..10}
    typer.run(enumerate_from_smiles)
