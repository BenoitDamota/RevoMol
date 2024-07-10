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


def list_neighbors(molecule: Molecule) -> set[Molecule]:
    molecule.compute_possible_actions()
    return {
        action.apply()
        for representation in molecule.possible_actions.values()
        for action_list in representation.values()
        for action in action_list
    }


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
) -> tuple[set[str], set[str]]:
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
    found_valid: set[str] = set()
    found_invalid: set[str] = set()

    if is_valid_molecule(Molecule(starting_smiles), evaluations):
        found_valid.add(starting_smiles)
    else:
        found_invalid.add(starting_smiles)

    while to_explore:
        # get the next SMILES to explore
        current_smiles: str = to_explore.pop()

        # find the neighbors of the current SMILES
        new_mols = find_neighbors(Molecule(current_smiles), depth)

        # check if the neighbors are valid and add them to the set to_try
        for new_mol in new_mols:
            new_smi = str(new_mol)
            if new_smi in found_valid or new_smi in found_invalid:
                continue

            if is_valid_molecule(new_mol, evaluations):
                to_explore.add(new_smi)
                found_valid.add(new_smi)
            else:
                found_invalid.add(new_smi)

    return found_valid, found_invalid


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

    print("molecule,heavy_atom_limit,depth,eval,nb_valid,nb_invalid,time")

    start_time = time.time()

    can_smi = Molecule(start_smiles).get_representation(MolecularGraph).canonical_smiles

    found_valid, found_invalid = explore_smiles(can_smi, depth, evaluations)

    duration = time.time() - start_time

    print(
        f"{can_smi},{nb_heavy_atoms},{depth},{eval_name},{len(found_valid)},{len(found_invalid)},{duration:.2f}"
    )

    file_path = f"output/enumeration_from_{can_smi}_{eval_name}_{nb_heavy_atoms}.txt"
    with open(file_path, "w") as file:
        for mol in found_valid:
            if str(mol):
                file.write(f"{mol}\n")
            else:
                file.write(f'""\n')


if __name__ == "__main__":
    # parallel -j 20 python enumerate_smiles.py C {1} {2} 2 ::: chembl chembl_zinc ::: {1..10}
    typer.run(enumerate_from_smiles)

    # enumerate_from_smiles("C", "chembl_zinc", 3, 2)
