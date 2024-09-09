"""
This script is used to find the closest valid neighbors of a molecule.
The molecule is represented by its SMILES string.
The validity of the neighbors is determined by the evaluation function
(use chembl or chembl_zinc).
You can also specify the max number of heavy atoms allowed.


examples from CLI:
python scripts/find_closest_neighbor.py C chembl 10

TTF
python scripts/find_closest_neighbor.py "C1=CSC(=C2SC=CS2)S1" chembl 10

DD
python scripts/find_closest_neighbor.py "N1=S=NC2=C1N=S=N2" chembl 10
"""

import os
import sys
import time

import typer

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol import evaluation as evaluator
from evomol.representation import MolecularGraph, Molecule
from evomol.search import enumeration as en


def find_closest_neighbors(
    start_smiles: str, filter_name: str, max_heavy_atoms: int
) -> None:
    """Look for the closest neighbors of a molecule that are valid and print them.

    Args:
        start_smiles (str): SMILES to start from
        eval_name (str): chembl or chembl_zinc as filter
        max_heavy_atoms (int): maximum number of heavy atoms
    """

    # init the default parameters and evaluations
    dp.setup_default_parameters(max_heavy_atoms=max_heavy_atoms)
    dp.setup_default_action_space()

    evaluations = dp.setup_filters(filter_name)

    # convert the starting SMILES to its canonical form
    can_smi_start = (
        Molecule(start_smiles).get_representation(MolecularGraph).canonical_smiles
    )

    print(f"{start_smiles} - canonical : {can_smi_start} - filter: {filter_name}")

    valid_mols: set[str] = set()
    depth: int = 0

    # explore the neighbors of the starting SMILES until a valid molecule is found
    while not valid_mols:
        depth += 1

        time_start = time.time()

        smiles_set: set[str] = en.find_neighbors(Molecule(can_smi_start), depth)

        for smi in smiles_set:
            if smi == can_smi_start:
                continue
            if smi in valid_mols:
                continue
            if evaluator.is_valid_molecule(Molecule(smi), evaluations):
                valid_mols.add(smi)

        duration = time.time() - time_start

        print(
            f"\t {len(valid_mols)}/{len(smiles_set)} valid molecules at "
            f"depth {depth} "
            f"in {duration:.2f}s"
        )

    for mol in valid_mols:
        print("\t\t", mol)


def main() -> None:
    """Search for the closest neighbors of a molecule that are valid."""
    smiles = [
        # "C1=CSC(=C2SC=CS2)S1",  # TTF
        # "N1=S=NC2=C1N=S=N2",  # DD
        # "CC1=NOC(C)(O)C1=NO",
        # "CN1ON(C)ON(C)O1",
        # "c1coc2occoc=2o1",
        # "CSC1N=NC(SC)N=N1",
        # "S=c1[nH]ssc2nnc1=2",
        # "N1=NC(=C2N=NN=N2)N=N1",
        # "N1=S=NC2=C1N=S=N2",
        "N#CN=C1SCS1",
        "C=CC1=CSNN1",
        "CSC1N=NN=N1",
    ]
    max_heavy_atoms = 10

    for smi in smiles:
        for eval_name in [
            "chembl",
            "chembl_zinc",
        ]:
            find_closest_neighbors(
                start_smiles=smi,
                filter_name=eval_name,
                max_heavy_atoms=max_heavy_atoms,
            )
        print()


if __name__ == "__main__":
    # main()

    typer.run(find_closest_neighbors)
