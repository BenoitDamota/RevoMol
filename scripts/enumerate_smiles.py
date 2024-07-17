"""
Short script to enumerate SMILES strings for a given molecule.

Parameters are:
- start_smiles: str - The SMILES string to start the enumeration from.
- eval_name: str - The name of the evaluation to use (chembl or chembl_zinc).
- nb_heavy_atoms: int - The max number of heavy atoms allowed.
- depth: int - The depth of the enumeration.

simple example:
python enumerate_smiles.py C chembl_zinc 3 2

parallel example:
parallel -j 20 python enumerate_smiles.py C {1} {2} 2 ::: chembl chembl_zinc ::: {1..10}
"""

import sys
import os

import typer

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


from evomol.search.enumeration import setup_and_launch_enumeration


if __name__ == "__main__":
    typer.run(setup_and_launch_enumeration)

    # setup_and_launch_enumeration("C", "chembl_zinc", 3, 2)
