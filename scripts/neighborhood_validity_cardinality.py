"""
This script is used to evaluate the cardinality of the neighborhood of a
molecule and the number of valid molecules in this neighborhood.
The script is used to compare the performance of the two functions :
- `find_neighbors_and_filter_with_duplicates`
- `find_neighbors_and_filter_without_duplicates`
from the evomol.search.enumeration module.
The script is used to compare the performance of the two functions in terms of 
the number of molecules in the neighborhood, the number of valid molecules in
this neighborhood and the time taken to compute the neighborhood and to filter 
the valid molecules in this neighborhood.

example:
    python scripts/neighborhood_validity_cardinality.py "CS(=O)(=O)C1=C(O)C(CC2=CC=C(F)S2)=CC(F)=C1"

gives the following output:

molecule,filter_between_depth,keep_duplicates,depth,nb_molecules,|molecules|,nb_valids,|valid|,time_generation,time_filter
CS(=O)(=O)C1=C(O)C(CC2=CC=C(F)S2)=CC(F)=C1,False,True,1,236,220,22,20,0.17,0.46
CS(=O)(=O)C1=C(O)C(CC2=CC=C(F)S2)=CC(F)=C1,False,True,2,58394,26708,1228,413,30.09,95.90
CS(=O)(=O)C1=C(O)C(CC2=CC=C(F)S2)=CC(F)=C1,False,False,1,220,220,20,20,0.13,0.35
CS(=O)(=O)C1=C(O)C(CC2=CC=C(F)S2)=CC(F)=C1,False,False,2,26708,26708,413,413,28.19,44.52

It means that starting from the molecule 
"CS(=O)(=O)C1=C(O)C(CC2=CC=C(F)S2)=CC(F)=C1" (which as "good results" in QED).
At a depth of 2 when not filtering molecules between each depth, if you keep
duplicate molecules in the neighborhood, you will have 58394 molecules to
explore only 26708 are unique (~= 46%) and only 413 are valid (~= 0.7%).
But if you ignore duplicate molecules in the neighborhood, you will have 26708
unique molecules to explore. Among these molecules, 413 are valid (~= 1.5%).

If apply evaluation is set to True, too many molecules are not reachable as
they are filtered between each depth.
"""

import os
import sys
import time
from typing import Collection

import typer

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol import evaluation as evaluator
from evomol.representation import Molecule
from evomol.search.enumeration import (
    find_neighbors_and_filter_with_duplicates,
    find_neighbors_and_filter_without_duplicates,
)


def neighborhood_validity_cardinality(smiles: str) -> None:
    """Print the cardinality of the neighborhood of a molecule with and without
    duplicates and the number of valid molecules in this neighborhood.

    Args:
        smiles (str): SMILES representation of the molecule.
    """
    dp.setup_default_parameters(
        accepted_atoms=["C", "O", "N", "F", "S"],
        # accepted_atoms=["C", "O", "N", "F", "S", "P", "Cl", "Br"],
        max_heavy_atoms=38,
    )
    dp.setup_default_action_space(
        with_add_group=False,
        with_remove_group=False,
    )

    evaluations = [
        evaluator.UnknownGCF(),
        evaluator.FilterUnknownGCF(threshold=0),
        evaluator.UnknownECFP(),
        evaluator.FilterUnknownECFP(threshold=0),
    ]

    print(
        "molecule,filter_between_depth,keep_duplicates,depth,"
        "nb_molecules,|molecules|,nb_valids,|valid|,time_generation,time_filter"
    )

    for apply_evaluation in [False]:
        for duplicates in [True, False]:
            for depth in range(1, 3):
                # explore the neighborhood of the molecule
                start = time.time()
                smiles_list: Collection[str] = (
                    find_neighbors_and_filter_with_duplicates(
                        Molecule(smiles), depth, evaluations, apply_evaluation
                    )
                    if duplicates
                    else find_neighbors_and_filter_without_duplicates(
                        Molecule(smiles), depth, evaluations, apply_evaluation
                    )
                )
                duration = time.time() - start

                # filter the valid molecules
                valid_mols = []
                start = time.time()
                for smi in smiles_list:
                    if evaluator.is_valid_molecule(Molecule(smi), evaluations):
                        valid_mols.append(smi)
                duration_post_filter = time.time() - start

                print(
                    ",".join(
                        map(
                            str,
                            [
                                smiles,
                                apply_evaluation,
                                duplicates,
                                depth,
                                len(smiles_list),
                                len(set(smiles_list)),
                                len(valid_mols),
                                len(set(valid_mols)),
                                f"{duration:.2f}",
                                f"{duration_post_filter:.2f}",
                            ],
                        )
                    )
                )
                # print(set(valid_mols))


if __name__ == "__main__":

    # smiles = "C"
    # smiles = "CCCCCCC"
    # smiles = "CS(=O)(=O)C1=C(O)C(CC2=CC=C(F)S2)=CC(F)=C1"
    # smiles = "C1=CSC(=C2SC=CS2)S1"  # TTF
    # smiles = "N1=S=NC2=C1N=S=N2"  # DD
    # smiles = (
    #     "C1=CC(=CC=C1C2=C3C=CC(=C(C4=NC(=C(C5=CC=C(N5)C(=C6C=CC2=N6)C7=CC="
    #     "C(C=C7)C(=O)O)C8=CC=C(C=C8)C(=O)O)C=C4)C9=CC=C(C=C9)C(=O)O)N3)C(=O)O"
    # )  # porph

    # neighborhood_validity_cardinality(smiles)

    typer.run(neighborhood_validity_cardinality)
