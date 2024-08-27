"""
This script launches an evolutionary algorithm to optimize the QED of a molecule.
"""

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol import evaluation as evaluator
from evomol.representation import MolecularGraph, Molecule
from evomol.search.evolutionary_algorithm import (
    EvolutionaryAlgorithm,
    ParametersEvolutionaryAlgorithm,
)
from evomol.search.neighborhood_strategy import RandomNeighborhoodStrategy
from evomol.search.parameters import search_parameters


def main() -> None:
    """Test the evolutionary algorithm."""

    # setup the default parameters
    dp.setup_default_parameters(
        accepted_atoms=["C", "O", "N", "F", "S"],
        max_heavy_atoms=38,
    )
    dp.setup_default_action_space(
        with_add_group=False,
        with_remove_group=False,
    )

    # the fitness criteria to optimize will be the QED
    search_parameters.fitness_criteria = "QED"

    # the evaluations to use
    evaluations = [
        evaluator.UnknownGCF(),
        evaluator.FilterUnknownGCF(threshold=0),
        evaluator.UnknownECFP(),
        evaluator.FilterUnknownECFP(threshold=0),
        evaluator.QED,
    ]

    # the population is composed of one empty molecule
    molecules = [Molecule("")]

    # ensure that the population is only canonical smiles
    population = [
        Molecule(molecule.get_representation(MolecularGraph).canonical_smiles)
        for molecule in molecules
    ]

    # evaluate the population
    for molecule in population:
        for evaluation in evaluations:
            evaluation.evaluate(molecule)

    # create the evolutionary algorithm
    evo_mol = EvolutionaryAlgorithm(
        population=population,
        parameters=ParametersEvolutionaryAlgorithm(),
        evaluations=evaluations,
        neighborhood_strategy=RandomNeighborhoodStrategy(
            depth=2,
            nb_max_tries=50,
            preselect_actions=True,
        ),
    )

    # launch the search
    evo_mol.run()


if __name__ == "__main__":
    main()
