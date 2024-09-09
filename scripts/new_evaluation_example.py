"""
Script to present the example of evaluation from the documentation.

The following example shows how to create an evaluation for a molecule with
three different cases :
1. The first case is when the evaluation is a function without any parameter or
    data (count atom and count bond).
2. The second case is when the evaluation is a class without any parameter or
    data but require other informations on the molecule (mean bonds per atom).
3. The second case is when the evaluation is a class with a parameter and data
    (filter mean bonds per atom).

"""

import math
import os
import sys

from typing_extensions import override

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# pylint: disable=wrong-import-position, import-error

import evomol.evaluation as evaluator
from evomol import default_parameters as dp
from evomol.evaluation import Evaluation, EvaluationError, Function, is_valid_molecule
from evomol.representation import MolecularGraph, Molecule

# 1. we define the functions to count the number of atoms and bonds in a
# molecule.
# They both use the MolecularGraph representation of the molecule to count the
# atoms and bonds.


def count_atoms(molecule: Molecule) -> int:
    """Count the number of atoms in the molecule."""
    return len(molecule.get_representation(MolecularGraph).atoms)


def count_bonds(molecule: Molecule) -> int:
    """Count the number of bonds in the molecule."""
    mol_graph = molecule.get_representation(MolecularGraph)
    return sum(mol_graph.bond_order(a1, a2) for a1, a2 in mol_graph.bonds)


# Then, we define the Evaluation objects for these functions.
CountAtoms = Function("count_atoms", count_atoms)
CountBonds = Function("count_bonds", count_bonds)
# "count_atoms" and "count_bonds" are the names of the evaluations, they are used
# to retrieve the values from the molecule object with the value method.

# 2. We define a function to evaluate the mean number of bonds per atom in a
# molecule. This function uses the count_atoms and count_bonds functions to
# compute the mean number of bonds per atom.


def mean_bonds_per_atom(molecule: Molecule) -> float:
    """Compute the mean number of bonds per atom in the molecule."""
    nb_atoms: int = molecule.value("count_atoms")
    nb_bonds: int = molecule.value("count_bonds")
    return nb_bonds / nb_atoms


# We define the Evaluation object for this function.
MeanBondsPerAtom = Function("mean_bonds_per_atom", mean_bonds_per_atom)


# 3. We define a class to filter molecules based on the mean number of bonds per
# atom. This class requires a threshold parameter to filter the molecules.
# To keep this value in memory, we create the FilterMeanBondsPerAtom class
# which inherits from Evaluation and implements the _evaluate method.
# It also gives a name to the evaluation, through the call to the super
# constructor.
# Other evaluations could need to load data from a file or a database, or to
# perform a complex computation that is not possible in a single function.
# This call also uses the EvaluationError exception to raise an error if the
# molecule does not pass the evaluation.


# pylint: disable=too-few-public-methods
class FilterMeanBondsPerAtom(Evaluation):
    """Filter molecules based on the mean number of bonds per atom."""

    def __init__(self, threshold: float):
        super().__init__("FilterMeanBondsPerAtom")
        self.threshold = threshold

    @override
    def _evaluate(self, molecule: Molecule) -> bool:
        value = molecule.value("mean_bonds_per_atom")

        if value <= self.threshold:
            return True

        raise EvaluationError(
            "The molecule has too many bonds per atom " f"({value} > {self.threshold})."
        )


# The main function creates a molecule from a SMILES string and evaluates it
# with the evaluations defined above.
# It prints the values of the evaluations and whether the molecule is valid or
# not.


def gaussian(x: float, mu: float, sig: float) -> float:
    """Compute the gaussian function."""
    return (
        1.0
        / (math.sqrt(2.0 * math.pi) * sig)
        * math.exp(-math.pow((x - mu) / sig, 2.0) / 2)
    )


def sigmoid(x: float) -> float:
    """Compute the sigmoid function."""
    return 1 / (1 + math.exp(-x))


def my_fitness(molecule: Molecule) -> float:
    """Compute the fitness of a molecule."""
    plogp: float = molecule.value("PLogP")
    qed: float = molecule.value("QED")
    value: float = gaussian(plogp, 10, 2) * sigmoid(qed)
    return value


def main() -> None:
    """Main function to evaluate molecules with the defined evaluations."""

    # first, we setup the default parameters (mainly for the representation)
    dp.setup_default_parameters()

    # We define the evaluations to use.
    # Note that the order of the evaluations is important, as the evaluations
    # are performed in the order they are given.
    # FilterMeanBondsPerAtom is the last evaluation, as it requires the
    # mean_bonds_per_atom value to be computed which needs the count_atoms and
    # count_bonds values to be computed.
    evaluations = [
        CountAtoms,  # need for MeanBondsPerAtom
        CountBonds,  # need for MeanBondsPerAtom
        MeanBondsPerAtom,  # need for FilterMeanBondsPerAtom
        FilterMeanBondsPerAtom(1),  # filter
        evaluator.QED,  # need for my_fitness
        evaluator.LogP,  # need for ZincNormalizedLogP
        evaluator.ZincNormalizedLogP,  # need for PLogP
        evaluator.SAScore,  # need for ZincNormalizedSAScore
        evaluator.ZincNormalizedSAScore,  # need for PLogP
        evaluator.CycleScore,  # need for NormalizedCycleScore
        evaluator.NormalizedCycleScore,  # need for PLogP
        evaluator.PLogP,  # need for my_fitness
        Function("my_fitness", my_fitness),
    ]

    # We define some SMILES to test the evaluations.
    smiles = [
        "CCC",
        "C=CC",
        "C=C=C",
        "C1=C=C=1",
    ]

    for smi in smiles:
        print(f"{smi:10} : ", end="")
        molecule = Molecule(smi)

        # We evaluate the molecule with the defined evaluations.
        result = is_valid_molecule(molecule, evaluations)
        if result:
            print("  valid")
        else:
            print("invalid")

        print(molecule.values)


# Output :
# CCC        :   valid
# {'count_atoms': 3, 'count_bonds': 2,
# 'mean_bonds_per_atom': 0.6666666666666666,
# 'FilterMeanBondsPerAtom': True, 'QED': 0.3854706587740357, 'logP': 1.4163,
# 'zinc_normalized_logP': -0.7256345487897407, 'SAScore': 1.7549570386349824,
# 'zinc_normalized_sa_score': 1.5567988735797864, 'CycleScore': 0,
# 'zinc_normalized_cycle_score': 0.16981148868758966,
# 'PLogP': 1.0009758134776354, 'my_fitness': 4.767151064685479e-06}
# C=CC       :   valid
# {'count_atoms': 3, 'count_bonds': 3, 'mean_bonds_per_atom': 1.0,
# 'FilterMeanBondsPerAtom': True, 'QED': 0.3570688167166915,
# 'logP': 1.1923, 'zinc_normalized_logP': -0.8818056352094792,
# 'SAScore': 3.006787807865753,
# 'zinc_normalized_sa_score': 0.054939632952350555,
# 'CycleScore': 0, 'zinc_normalized_cycle_score': 0.16981148868758966,
# 'PLogP': -0.6570545135695389, 'my_fitness': 8.016705725675326e-08}
# C=C=C      : invalid
# {'count_atoms': 3, 'count_bonds': 4,
# 'mean_bonds_per_atom': 1.3333333333333333}
# C1=C=C=1   : invalid
# {'count_atoms': 3, 'count_bonds': 6, 'mean_bonds_per_atom': 2.0}

if __name__ == "__main__":
    main()
