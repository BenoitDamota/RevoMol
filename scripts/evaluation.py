"""
This script demonstrates how to use the evaluation module to evaluate molecules.
It defines a custom fitness function and uses it to evaluate molecules.

error_order function demonstrates an error when the order of evaluations is
incorrect, the custom fitness function is evaluated before logP and QED.
"""

import os
import sys
from pprint import pprint

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error

import evomol.default_parameters as dp
from evomol import evaluation as evaluator
from evomol.representation import Molecule


def my_fitness(molecule: Molecule) -> float:
    """Custom fitness function that combines logP and QED.

    Args:
        molecule (Molecule): Molecule to evaluate.

    Returns:
        float: Fitness value.
    """
    log_p: float = molecule.value("logP")
    qed: float = molecule.value("QED")
    return log_p + qed


def main() -> None:
    """Example of using the evaluation module to evaluate molecules."""
    print("Good evaluation example")
    dp.setup_default_parameters()

    # Define evaluations
    # UnknownGCF and UnknownECFP are used to check if the molecule is valid
    # FilterUnknownGCF and FilterUnknownECFP are used to filter out if
    # the molecule has unknown GCF or ECFP
    # LogP and QED are used to evaluate the molecule
    # my_fitness is a custom fitness function that combines logP and QED
    # and has the name "my_fitness"
    evaluations = [
        evaluator.UnknownGCF(),
        evaluator.FilterUnknownGCF(),
        evaluator.UnknownECFP(),
        evaluator.FilterUnknownECFP(),
        evaluator.LogP,
        evaluator.QED,
        evaluator.Function("my_fitness", my_fitness),
    ]

    # SMILES to evaluate
    smiles = [
        "C",
        "C=1=CC1",  # this molecule is invalid (unknown GCF)
        "c1ccccc1",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(=O)NC1=CC=C(C=C1)O",
        "C1C2C1CC2",
    ]

    for smi in smiles:
        mol = Molecule(smi)
        is_valid = evaluator.is_valid_molecule(mol, evaluations)
        print(f"{smi} is{'' if is_valid else ' not'} valid")
        pprint(mol.values)
        print()


def error_order() -> None:
    """
    Example of error when the order of evaluations is incorrect.
    If the custom fitness function is evaluated before logP and QED,
    an error KeyError will be raised because logP and QED are not evaluated yet.
    This is a critical error so the program will stop.
    """
    print("Error evaluation example")
    dp.setup_default_parameters()

    evaluations = [
        evaluator.UnknownGCF(),
        evaluator.FilterUnknownGCF(),
        evaluator.UnknownECFP(),
        evaluator.FilterUnknownECFP(),
        # error here because my_fitness is evaluated before logP and QED
        evaluator.Function("my_fitness", my_fitness),
        evaluator.LogP,
        evaluator.QED,
    ]

    smiles = [
        "C",
        "C=1=CC1",  # this molecule is invalid (unknown GCF)
        "c1ccccc1",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
    ]

    for smi in smiles:
        mol = Molecule(smi)
        is_valid = evaluator.is_valid_molecule(mol, evaluations)
        print(f"{smi} is{'' if is_valid else ' not'} valid")
        pprint(mol.values)
        print()


if __name__ == "__main__":
    main()
    print(f"\n{'='*50}\n")
    error_order()
