"""
Short script to enumerate SMILES strings for a given molecule.

Parameters are:
- start_smiles: str - The SMILES string to start the enumeration from.
- eval_name: str - The name of the evaluation to use (chembl or chembl_zinc).
- nb_heavy_atoms: int - The max number of heavy atoms allowed.
- depth: int - The depth of the enumeration.
- nb_cpu: int - The number of CPU to use. 0 to use all CPUs.
- enumeration_type: str - The type of enumeration to use (sequential, 
    multiprocessing or workers).

simple example:
python scripts/enumerate_smiles.py C chembl_zinc 3 2 0 multiprocessing

parallel example (join the 2 lines):
parallel -j 20 python scripts/enumerate_smiles.py 
    C {1} {2} 2 0 multiprocessing ::: chembl chembl_zinc ::: {1..10}
"""

import os
import sys
import time

import typer

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol.representation import MolecularGraph, Molecule
from evomol.search import enumeration as en


def setup_and_launch_enumeration(
    start_smiles: str,
    eval_name: str,
    nb_heavy_atoms: int,
    depth: int,
    nb_cpu: int,
    enumeration_type: str,
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
        nb_cpu (int): number of CPU to use. 0 to use all CPUs.
        enumeration_type (str): sequential, multiprocessing or workers.
            sequential uses a simple loop to explore the neighbors.
            multiprocessing uses the multiprocessing library to parallelize the
            exploration of the neighbors (not the filtering).
            workers uses a queue of smiles to explore that workers will take
            from to explore the neighbors (not the filtering).
            Defaults to parallel.
    """

    # init the default parameters and evaluations
    dp.setup_default_parameters(max_heavy_atoms=nb_heavy_atoms)
    dp.setup_default_action_space()
    evaluations = dp.setup_filters(eval_name)

    # print the header of the CSV
    print("molecule,heavy_atom_limit,depth,eval,nb_valid,nb_invalid,time")

    start_time = time.time()

    # convert the starting SMILES to its canonical form
    can_smi = Molecule(start_smiles).get_representation(MolecularGraph).canonical_smiles

    # explore the neighbors of the starting SMILES
    found_valid: set[str]
    found_invalid: set[str]

    if enumeration_type == "sequential":
        found_valid, found_invalid = en.enumerate_from_smiles(
            can_smi, depth, evaluations
        )

    if enumeration_type == "multiprocessing":
        found_valid, found_invalid = en.enumerate_from_smiles_multiprocessing(
            can_smi, depth, evaluations, nb_cpu
        )
    if enumeration_type == "workers":
        found_valid, found_invalid = en.enumerate_from_smiles_workers(
            can_smi, depth, evaluations, nb_cpu
        )

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
            str_mol: str = str(Molecule(mol))
            if str_mol:
                file.write(f"{str_mol}\n")
            # don't write the empty molecule
            # else:
            #     file.write('""\n')


if __name__ == "__main__":
    typer.run(setup_and_launch_enumeration)

    # setup_and_launch_enumeration("C", "chembl_zinc", 10, 2, 0, "multiprocessing")
    # setup_and_launch_enumeration("C", "chembl_zinc", 4, 2, 0, "workers")
