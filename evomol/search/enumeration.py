"""
Those functions are used to explore the neighbors of a molecule up to a certain depth.
Either by listing all the neighbors for one molecule and printing the valid ones
with list_up_to_depth_from_smiles or by exploring the complete neighborhood of 
a molecule and the neighborhood of its neighbors and so on with 
enumerate_from_smiles.

Those functions could be optimized.

find_neighbors could take the set of valid neighbors already found and only
consider the neighbors that are not in this set, I tried it but it was slower
and didn't lead to the same results. Maybe I did something wrong.

enumerate_from_smiles could be optimized by using multiprocessing to explore
the neighbors of the molecules in parallel.
Two pools of workers could be used, one to explore the neighbors of the molecules
and the other to check the validity of the molecules.
It tried to implement it but I had some issues with the Molecule class it seems
that it didn't keep the default values of the class attributes. I didn't have
time to investigate further.
"""

import time

from evomol import evaluation as evaluator
from evomol.action import molecular_graph as mg
from evomol.evaluation import is_valid_molecule
from evomol.representation import SMILES, MolecularGraph, Molecule


def list_and_compute_neighbors(molecule: Molecule) -> set[Molecule]:
    """List the neighbors of a molecule and compute them

    Args:
        molecule (Molecule): Molecule to explore

    Returns:
        set[Molecule]: set of molecules found
    """
    # fill the possible actions of the molecule
    molecule.compute_possible_actions()
    # for each representation, for each action list, for each action, apply it
    # and return the set of new molecules
    return {
        action.apply()
        for representation in molecule.possible_actions.values()
        for action_list in representation.values()
        for action in action_list
    }


def find_neighbors(molecule: Molecule, max_depth: int) -> set[Molecule]:
    """Iteratively find the neighbors of a molecule up to a certain depth.
    Explore all neighbors at level 1 then all neighbors at level 2 and so on.

    Args:
        molecule (Molecule): Molecule to explore
        max_depth (int): Maximum depth to explore

    Returns:
        set[Molecule]: set of molecules found
    """
    # set of neighbors found
    neighbors: set[Molecule] = set()

    # queue of molecules to explore
    queue = {molecule}

    # start at depth 1
    depth = 1

    # explore the neighbors of the molecule up to the maximum depth
    while depth <= max_depth:
        # set of neighbors found in this depth
        next_queue = set()
        # for each molecule in the queue, find the neighbors
        for current_mol in queue:
            for new_mol in list_and_compute_neighbors(current_mol):
                # add the new molecule to the set of neighbors if it is not
                # already in it, don't add it if it is already in the set
                # as it has already been explored
                if new_mol not in neighbors:
                    next_queue.add(new_mol)
                    neighbors.add(new_mol)
        # update the queue with the new neighbors and increase the depth
        queue = next_queue
        depth += 1

    return neighbors


def find_neighbors_and_check_validity_each_step(
    molecule: Molecule, max_depth: int, evaluations: list[evaluator.Evaluation]
) -> set[Molecule]:
    """Iteratively find the neighbors of a molecule up to a certain depth while
    checking the validity of each molecule at each step.
    Same as find_neighbors but check the validity of the molecules at each step.

    Args:
        molecule (Molecule): Molecule to explore
        max_depth (int): Maximum depth to explore

    Returns:
        set[Molecule]: set of molecules found
    """
    neighbors: set[Molecule] = set()
    queue: set[Molecule] = {molecule}
    depth: int = 1

    while depth <= max_depth:
        next_queue: set[Molecule] = set()
        for current_mol in queue:
            for new_mol in list_and_compute_neighbors(current_mol):
                if new_mol not in neighbors and is_valid_molecule(new_mol, evaluations):
                    next_queue.add(new_mol)
                    neighbors.add(new_mol)
        queue = next_queue
        depth += 1

    return neighbors


def enumerate_from_smiles(
    starting_smiles: str, depth: int, evaluations: list[evaluator.Evaluation]
) -> tuple[set[str], set[str]]:
    """Enumerate all reachable molecules starting from a SMILES up to a certain
    depth and check their validity before exploring them.

    Args:
        starting_smiles (str): SMILES to explore
        depth (int): Maximum depth to explore
        evaluations (list[evaluator.Evaluation]): filter to apply

    Returns:
        tuple[set[str], set[str]]: sets of valid and invalid SMILES
    """

    # set of SMILES to explore
    to_explore: set[str] = {starting_smiles}

    # set of valid and invalid SMILES found
    found_valid: set[str] = set()
    found_invalid: set[str] = set()

    # add the starting SMILES to the set of valid or invalid SMILES
    if is_valid_molecule(Molecule(starting_smiles), evaluations):
        found_valid.add(starting_smiles)
    else:
        found_invalid.add(starting_smiles)

    while to_explore:
        # get the next SMILES to explore
        current_smiles: str = to_explore.pop()

        # find the neighbors of the current SMILES
        new_mols = find_neighbors(Molecule(current_smiles), depth)

        for new_mol in new_mols:
            # get the SMILES of the new molecule
            new_smi = str(new_mol)
            # if the SMILES is already in the set of valid or invalid SMILES
            # it has already been explored or is in the queue to be explored
            if new_smi in found_valid or new_smi in found_invalid:
                continue

            # check if the new molecule is valid
            if is_valid_molecule(new_mol, evaluations):
                # if it is valid, add it to the set of valid SMILES and to the
                # set of SMILES to explore
                to_explore.add(new_smi)
                found_valid.add(new_smi)
            else:
                # if it is invalid, add it to the set of invalid SMILES
                # to avoid checking it again
                found_invalid.add(new_smi)

    return found_valid, found_invalid


def evaluation_gcf_ecfp(eval_name: str) -> list[evaluator.Evaluation]:
    """Return the evaluations for the generic cyclic features and ECFP

    Args:
        eval_name (str): chembl or chembl_zinc as filter

    Returns:
        list[evaluator.Evaluation]: list of evaluations
    """
    if eval_name == "chembl":
        return [
            evaluator.UnknownGCF(path_db="external_data/gcf1.txt", name="chembl"),
            evaluator.FilterUnknownGCF(threshold=0, name="chembl"),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL.txt", radius=2, name="chembl"
            ),
            evaluator.FilterUnknownECFP(threshold=0, name="chembl"),
        ]
    if eval_name == "chembl_zinc":
        return [
            evaluator.UnknownGCF(path_db="external_data/gcf2.txt", name="chembl_zinc"),
            evaluator.FilterUnknownGCF(threshold=0, name="chembl_zinc"),
            evaluator.UnknownECFP(
                path_db="external_data/ecfp4_ChEMBL_ZINC.txt",
                radius=2,
                name="chembl_zinc",
            ),
            evaluator.FilterUnknownECFP(threshold=0, name="chembl_zinc"),
        ]
    raise ValueError("eval_name must be 'chembl' or 'chembl_zinc'")


def setup_default_enumeration_parameters(nb_heavy_atoms: int) -> None:
    """Set the default parameters for the enumeration

    Args:
        nb_heavy_atoms (int): maximum number of heavy atoms
    """
    Molecule.id_representation_class = SMILES
    Molecule.representations_class = [MolecularGraph]
    Molecule.accepted_atoms = ["C", "O", "N", "F", "S"]

    Molecule.max_heavy_atoms = nb_heavy_atoms

    # only use the "classical" actions
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


def setup_and_launch_enumeration(
    start_smiles: str, eval_name: str, nb_heavy_atoms: int, depth: int
) -> None:
    """Setup the default parameters and launch the enumeration from a SMILES.
    All valid molecules found are saved in a file in the output/ directory,
    one molecule per line :
    enumeration_{can_smi}_{eval_name}_max_atom_{nb_heavy_atoms}_depth_{depth}.txt

    The time taken to explore the molecules and the number of valid and invalid
    molecules found is printed in a CSV format.

    Args:
        start_smiles (str): SMILES to start from
        eval_name (str): chembl or chembl_zinc as filter
        nb_heavy_atoms (int): maximum number of heavy atoms
        depth (int): maximum depth to explore
    """

    # init the default parameters and evaluations
    setup_default_enumeration_parameters(nb_heavy_atoms)
    evaluations = evaluation_gcf_ecfp(eval_name)

    # print the header of the CSV
    print("molecule,heavy_atom_limit,depth,eval,nb_valid,nb_invalid,time")

    start_time = time.time()

    # convert the starting SMILES to its canonical form
    can_smi = Molecule(start_smiles).get_representation(MolecularGraph).canonical_smiles

    # explore the neighbors of the starting SMILES
    found_valid, found_invalid = enumerate_from_smiles(can_smi, depth, evaluations)

    duration = time.time() - start_time

    # print the results in a CSV format
    print(
        f"{can_smi},{nb_heavy_atoms},{depth},{eval_name},"
        f"{len(found_valid)},{len(found_invalid)},{duration:.2f}"
    )

    # print the valid molecules in a file, one molecule per line
    file_path = (
        f"output/enumeration_{can_smi}_{eval_name}"
        f"_max_atom_{nb_heavy_atoms}_depth_{depth}.txt"
    )
    with open(file_path, "w", encoding="utf8") as file:
        for mol in found_valid:
            if str(mol):
                file.write(f"{mol}\n")
            else:
                file.write('""\n')


def list_up_to_depth_from_smiles(
    start_smiles: str, eval_name: str, nb_heavy_atoms: int, depth: int
) -> None:
    """List the neighbors of a molecule up to a certain depth and print the valid ones.

    Args:
        start_smiles (str): SMILES to start from
        eval_name (str): chembl or chembl_zinc as filter
        nb_heavy_atoms (int): maximum number of heavy atoms
        depth (int): maximum depth to explore
    """

    # init the default parameters and evaluations
    setup_default_enumeration_parameters(nb_heavy_atoms)
    evaluations = evaluation_gcf_ecfp(eval_name)

    # convert the starting SMILES to its canonical form
    can_smi = Molecule(start_smiles).get_representation(MolecularGraph).canonical_smiles

    print(f"Starting from {can_smi} (from {start_smiles})")

    # find the neighbors of the starting SMILES
    neighbors = find_neighbors(Molecule(can_smi), depth)

    # print the valid neighbors
    print("Neighbors:")
    for mol in neighbors:
        if is_valid_molecule(mol, evaluations):
            print("\t", mol)


def find_closest_neighbors(
    start_smiles: str, eval_name: str, nb_heavy_atoms: int
) -> None:
    """Look for the closest neighbors of a molecule that are valid and print them.

    Args:
        start_smiles (str): SMILES to start from
        eval_name (str): chembl or chembl_zinc as filter
        nb_heavy_atoms (int): maximum number of heavy atoms
    """

    # init the default parameters and evaluations
    setup_default_enumeration_parameters(nb_heavy_atoms)
    evaluations = evaluation_gcf_ecfp(eval_name)

    # convert the starting SMILES to its canonical form
    can_smi_start = (
        Molecule(start_smiles).get_representation(MolecularGraph).canonical_smiles
    )

    print(f"{start_smiles} - canonical : {can_smi_start} - filter: {eval_name}")

    valid_mols: set[str] = set()
    depth: int = 0

    # explore the neighbors of the starting SMILES until a valid molecule is found
    while not valid_mols:
        depth += 1

        time_start = time.time()

        mols: set[Molecule] = find_neighbors(Molecule(can_smi_start), depth)

        for mol in mols:
            can_smi = str(mol)
            if can_smi == can_smi_start:
                continue
            if can_smi in valid_mols:
                continue
            if is_valid_molecule(mol, evaluations):
                valid_mols.add(str(mol))

        duration = time.time() - time_start

        print(
            f"\t {len(valid_mols)}/{len(mols)} valid molecules at "
            f"depth {depth} "
            f"in {duration:.2f}s"
        )

    for mol in valid_mols:
        print("\t\t", mol)
