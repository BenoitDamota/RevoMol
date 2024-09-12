"""
Those functions are used to explore the neighbors of a molecule up to a certain depth.
Either by listing all the neighbors for one molecule and printing the valid ones
with list_up_to_depth_from_smiles or by exploring the complete neighborhood of
a molecule and the neighborhood of its neighbors and so on with
enumerate_from_smiles.

Those functions may be be optimized.

find_neighbors could take the set of valid neighbors already found and only
consider the neighbors that are not in this set, I tried it but it was slower
and didn't lead to the same results. Maybe I did something wrong.

enumerate_from_smiles_parallel is optimized by using multiprocessing to explore
the neighbors of the molecules in parallel.
I tried to use multiprocessing to check the validity of the molecules but it was
slower :
# to_check: set[str] = set()
# for result in neighbors_results:
#     for new_smi in result:
#         # if the SMILES is already in the set of valid or invalid SMILES
#         # it has already been explored or is in the queue to be explored
#         if new_smi in found_valid or new_smi in found_invalid:
#             continue

#         to_check.add(new_smi)

# validity_results = pool.starmap(
#     parallel_is_valid_molecule, [(smi, evaluations) for smi in to_check]
# )
# for new_smi, is_valid in validity_results:
#     if is_valid:
#         to_explore.add(new_smi)
#         found_valid.add(new_smi)
#     else:
#         found_invalid.add(new_smi)
"""

from multiprocessing import Manager, Pool, Queue, cpu_count
from multiprocessing.managers import ValueProxy
from threading import Lock

from evomol import evaluation as evaluator
from evomol.representation import Molecule


def set_of_neighbors(molecule: Molecule) -> set[Molecule]:
    """List the neighbors of a molecule and compute them.

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


def list_of_neighbors(molecule: Molecule) -> list[Molecule]:
    """List the neighbors of a molecule and compute them.

    Args:
        molecule (Molecule): Molecule to explore

    Returns:
        list[Molecule]: list of molecules found
    """
    # fill the possible actions of the molecule
    molecule.compute_possible_actions()
    # for each representation, for each action list, for each action, apply it
    # and return the list of new molecules
    return [
        action.apply()
        for representation in molecule.possible_actions.values()
        for action_list in representation.values()
        for action in action_list
    ]


def dict_of_neighbors(molecule: Molecule) -> dict[str, dict[str, list[Molecule]]]:
    """List the neighbors of a molecule and compute them.

    Args:
        molecule (Molecule): Molecule to explore

    Returns:
        dict[str, list[Molecule]]: dict of molecules found
    """
    # fill the possible actions of the molecule
    molecule.compute_possible_actions()
    # for each representation, for each action list, for each action, apply it
    # and return the dict of new molecules
    return {
        representation_name: {
            action_name: [action.apply() for action in action_list]
            for action_name, action_list in action_dict.items()
        }
        for representation_name, action_dict in molecule.possible_actions.items()
    }


def find_neighbors_and_filter_with_duplicates(
    molecule: Molecule,
    max_depth: int,
    evaluations: list[evaluator.Evaluation],
    apply_evaluation: bool,
) -> list[str]:
    """Iteratively find the neighbors of a molecule up to a certain depth.
    Explore all neighbors at level 1 then all neighbors at level 2 and so on.
    Duplicate molecules are kept.
    Molecules can be filtered with the evaluations between each depth.

    Args:
        molecule (Molecule): Molecule to explore
        max_depth (int): Maximum depth to explore

    Returns:
        list[str]: list of molecules found
    """
    # list of neighbors found
    neighbors: list[str] = []

    # queue of molecules to explore
    queue: list[Molecule] = [molecule]

    # start at depth 1
    depth = 1

    # explore the neighbors of the molecule up to the maximum depth
    while depth <= max_depth:
        # set of neighbors found in this depth
        next_queue: list[Molecule] = []
        # for each molecule in the queue, find the neighbors
        for current_mol in queue:
            for new_mol in list_of_neighbors(current_mol):
                new_smi = str(new_mol)
                if apply_evaluation and not evaluator.is_valid_molecule(
                    new_mol, evaluations
                ):
                    continue
                next_queue.append(new_mol)
                neighbors.append(new_smi)
        # update the queue with the new neighbors and increase the depth
        queue = next_queue
        depth += 1

    return neighbors


def find_neighbors_and_filter_without_duplicates(
    molecule: Molecule,
    max_depth: int,
    evaluations: list[evaluator.Evaluation],
    apply_evaluation: bool,
) -> set[str]:
    """Iteratively find the neighbors of a molecule up to a certain depth.
    Explore all neighbors at level 1 then all neighbors at level 2 and so on.
    Duplicate molecules are removed.
    Molecules can be filtered with the evaluations between each depth.

    Args:
        molecule (Molecule): Molecule to explore
        max_depth (int): Maximum depth to explore

    Returns:
        list[str]: list of molecules found
    """
    # set of neighbors found
    neighbors: set[str] = set()

    # queue of molecules to explore
    queue: set[Molecule] = {molecule}

    # start at depth 1
    depth = 1

    # explore the neighbors of the molecule up to the maximum depth
    while depth <= max_depth:
        # set of neighbors found in this depth
        next_queue: set[Molecule] = set()
        # for each molecule in the queue, find the neighbors
        for current_mol in queue:
            for new_mol in list_of_neighbors(current_mol):
                new_smi = str(new_mol)
                if apply_evaluation and not evaluator.is_valid_molecule(
                    new_mol, evaluations
                ):
                    continue
                if new_smi not in neighbors:
                    next_queue.add(new_mol)
                    neighbors.add(new_smi)
        # update the queue with the new neighbors and increase the depth
        queue = next_queue
        depth += 1

    return neighbors


def find_neighbors(molecule: Molecule, max_depth: int) -> set[str]:
    """Iteratively find the neighbors of a molecule up to a certain depth.
    Explore all neighbors at level 1 then all neighbors at level 2 and so on.

    Args:
        molecule (Molecule): Molecule to explore
        max_depth (int): Maximum depth to explore

    Returns:
        set[str]: set of molecules found
    """
    # set of neighbors found
    neighbors: set[str] = set()

    # queue of molecules to explore
    queue: set[Molecule] = {molecule}

    # start at depth 1
    depth = 1

    # explore the neighbors of the molecule up to the maximum depth
    while depth <= max_depth:
        # set of neighbors found in this depth
        next_queue = set()
        # for each molecule in the queue, find the neighbors
        for current_mol in queue:
            for new_mol in set_of_neighbors(current_mol):
                new_smi = str(new_mol)
                # add the new molecule to the set of neighbors if it is not
                # already in it, don't add it if it is already in the set
                # as it has already been explored
                if new_smi not in neighbors:
                    next_queue.add(new_mol)
                    neighbors.add(new_smi)
        # update the queue with the new neighbors and increase the depth
        queue = next_queue
        depth += 1

    return neighbors


def parallel_find_neighbors(smiles: str, depth: int) -> set[str]:
    """Find the neighbors of a molecule up to a certain depth for parallel use.

    Args:
        smiles (str): SMILES to explore
        depth (int): Maximum depth to explore

    Returns:
        set[str]: set of molecules found
    """
    return find_neighbors(Molecule(smiles), depth)


def parallel_is_valid_molecule(
    smiles: str, evaluations: list[evaluator.Evaluation]
) -> tuple[str, bool]:
    """Check if a molecule is valid according to a list of evaluations for
    parallel use.

    Args:
        smiles (str): SMILES to check
        evaluations (list[evaluator.Evaluation]): filter to apply

    Returns:
        tuple[str, bool]: _description_
    """
    return smiles, evaluator.is_valid_molecule(Molecule(smiles), evaluations)


def worker_find_neighbors(
    explore_queue: Queue,
    neighbor_queue: Queue,
    depth: int,
    task_counter: ValueProxy[int],
    lock: Lock,
) -> None:
    """Worker to explore the neighbors of a molecule up to a certain depth.

    Args:
        explore_queue (Queue): queue of SMILES to explore
        neighbor_queue (Queue): queue of neighbors found
        depth (int): Maximum depth to explore
        task_counter (Value): shared counter of tasks
    """
    while True:
        smiles = explore_queue.get()

        neighbors = find_neighbors(Molecule(smiles), depth)
        neighbor_queue.put(neighbors)
        with lock:
            task_counter.value -= 1


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
    if evaluator.is_valid_molecule(Molecule(starting_smiles), evaluations):
        found_valid.add(starting_smiles)
    else:
        found_invalid.add(starting_smiles)

    while to_explore:

        to_check: set[str] = set()

        while to_explore:
            # get the next SMILES to explore
            current_smiles: str = to_explore.pop()

            # find the neighbors of the current SMILES
            new_mols = find_neighbors(Molecule(current_smiles), depth)

            to_check.update(new_mols)

        while to_check:
            new_smi = to_check.pop()

            # if the SMILES is already in the set of valid or invalid SMILES
            # it has already been explored or is in the queue to be explored
            if new_smi in found_valid or new_smi in found_invalid:
                continue

            # check if the new molecule is valid
            if evaluator.is_valid_molecule(Molecule(new_smi), evaluations):
                # if it is valid, add it to the set of valid SMILES and to the
                # set of SMILES to explore
                to_explore.add(new_smi)
                found_valid.add(new_smi)
            else:
                # if it is invalid, add it to the set of invalid SMILES
                # to avoid checking it again
                found_invalid.add(new_smi)

    return found_valid, found_invalid


def enumerate_from_smiles_multiprocessing(
    starting_smiles: str,
    depth: int,
    evaluations: list[evaluator.Evaluation],
    nb_cpu: int,
) -> tuple[set[str], set[str]]:
    """Enumerate all reachable molecules starting from a SMILES up to a certain
    depth and check their validity before exploring them in parallel.

    Args:
        starting_smiles (str): SMILES to explore
        depth (int): Maximum depth to explore
        evaluations (list[evaluator.Evaluation]): filter to apply
        nb_cpu (int): number of CPU to use

    Returns:
        tuple[set[str], set[str]]: sets of valid and invalid SMILES
    """

    # set of SMILES to explore
    to_explore: set[str] = {starting_smiles}

    # set of valid and invalid SMILES found
    found_valid: set[str] = set()
    found_invalid: set[str] = set()

    # add the starting SMILES to the set of valid or invalid SMILES
    if evaluator.is_valid_molecule(Molecule(starting_smiles), evaluations):
        found_valid.add(starting_smiles)
    else:
        found_invalid.add(starting_smiles)

    if nb_cpu == 0:
        nb_cpu = cpu_count()

    with Pool(nb_cpu) as pool:
        while to_explore:

            # find the neighbors of the SMILES to explore in parallel
            neighbors_results = pool.starmap(
                parallel_find_neighbors,
                [(smiles, depth) for smiles in to_explore],
            )
            to_explore.clear()

            for result in neighbors_results:
                for new_smi in result:
                    # if the SMILES is already in the set of valid or invalid SMILES
                    # it has already been explored or is in the queue to be explored
                    if new_smi in found_valid or new_smi in found_invalid:
                        continue

                    # check if the new molecule is valid
                    if evaluator.is_valid_molecule(Molecule(new_smi), evaluations):
                        # if it is valid, add it to the set of valid SMILES and to the
                        # set of SMILES to explore
                        to_explore.add(new_smi)
                        found_valid.add(new_smi)
                    else:
                        # if it is invalid, add it to the set of invalid SMILES
                        # to avoid checking it again
                        found_invalid.add(new_smi)

    return found_valid, found_invalid


def enumerate_from_smiles_workers(
    starting_smiles: str,
    depth: int,
    evaluations: list[evaluator.Evaluation],
    nb_cpu: int,
) -> tuple[set[str], set[str]]:
    """Enumerate all reachable molecules starting from a SMILES up to a certain
    depth and check their validity before exploring them in parallel.

    Args:
        starting_smiles (str): SMILES to explore
        depth (int): Maximum depth to explore
        evaluations (list[evaluator.Evaluation]): filter to apply
        nb_cpu (int): number of CPU to use

    Returns:
        tuple[set[str], set[str]]: sets of valid and invalid SMILES
    """

    if nb_cpu == 0:
        nb_cpu = cpu_count()

    manager = Manager()
    explore_queue = manager.Queue()
    neighbor_queue = manager.Queue()
    lock = manager.Lock()
    task_counter = manager.Value(int, 0)
    found_valid = set()
    found_invalid = set()

    # add the starting SMILES to the set of valid or invalid SMILES
    if evaluator.is_valid_molecule(Molecule(starting_smiles), evaluations):
        found_valid.add(starting_smiles)
    else:
        found_invalid.add(starting_smiles)

    with lock:
        task_counter.value += 1
    # starting SMILES is the first molecule to explore
    explore_queue.put(starting_smiles)

    # start the workers
    with Pool(nb_cpu) as pool:
        for _ in range(nb_cpu):
            pool.apply_async(
                worker_find_neighbors,
                (explore_queue, neighbor_queue, depth, task_counter, lock),
            )

        while True:
            while not neighbor_queue.empty():
                new_smis = neighbor_queue.get()
                for new_smi in new_smis:
                    if new_smi in found_valid or new_smi in found_invalid:
                        continue

                    if evaluator.is_valid_molecule(Molecule(new_smi), evaluations):
                        explore_queue.put(new_smi)
                        with lock:
                            task_counter.value += 1
                        found_valid.add(new_smi)
                    else:
                        found_invalid.add(new_smi)

            # if all tasks are done, break the loop
            if task_counter.value == 0:
                break

    return found_valid, found_invalid
