"""
This script is used to find the closest valid neighbors of a molecule.
The molecule is represented by its SMILES string.
The validity of the neighbors is determined by the evaluation function
(use chembl or chembl_zinc).
You can also specify the max number of heavy atoms allowed.


examples from CLI:
python find_closest_neighbor.py C chembl 10

TTF
python find_closest_neighbor.py "C1=CSC(=C2SC=CS2)S1" chembl 10

DD
python find_closest_neighbor.py "N1=S=NC2=C1N=S=N2" chembl 10
"""

import os
import sys

import typer

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from evomol.search.enumeration import find_closest_neighbors


def main():

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
    nb_heavy_atoms = 10

    for smi in smiles:
        for eval_name in [
            "chembl",
            "chembl_zinc",
        ]:
            find_closest_neighbors(
                start_smiles=smi,
                eval_name=eval_name,
                nb_heavy_atoms=nb_heavy_atoms,
            )
        print()


if __name__ == "__main__":
    # main()

    typer.run(find_closest_neighbors)
