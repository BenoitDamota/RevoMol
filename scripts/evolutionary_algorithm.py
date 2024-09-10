"""
This script launches an evolutionary algorithm to maximize the QED of a molecule.
"""

import os
import random
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error

from evomol import default_parameters as dp
from evomol import evaluation as evaluator
from evomol.logging import init_logger
from evomol.representation import MolecularGraph, Molecule
from evomol.search import selection_strategy
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

    logger = init_logger(
        config_file=os.path.join("logging_configs", "stderr-json-file.json"),
        output_file=os.path.join("logs", "evomol_default.log"),
    )

    random_seed = random.randint(0, 1000)
    logger.info(
        "Init the evolutionary algorithm.",
        extra={"random_seed": random_seed},
    )

    random.seed(random_seed)

    # the fitness criteria to optimize will be the QED
    search_parameters.fitness_criteria = "QED"

    # the evaluations to use
    evaluations = [
        evaluator.UnknownGCF(
            path_db=os.path.join("external_data", "gcf1.txt"),
            name="chembl",
        ),
        evaluator.FilterUnknownGCF(threshold=0, name="chembl"),
        evaluator.UnknownECFP(
            path_db=os.path.join("external_data", "ecfp4_ChEMBL.txt"),
            radius=2,
            name="chembl",
        ),
        evaluator.FilterUnknownECFP(threshold=0, name="chembl"),
        evaluator.QED,
    ]

    # the population is composed of one empty molecule
    molecules = [Molecule("C")]

    # ensure that the population is only canonical smiles
    population = [
        Molecule(molecule.get_representation(MolecularGraph).canonical_smiles)
        for molecule in molecules
    ]

    # create the evolutionary algorithm
    evo_mol = EvolutionaryAlgorithm(
        population=population,
        parameters=ParametersEvolutionaryAlgorithm(
            size=1000,
            max_generations=1000,
            max_time=3600,
            nb_replacements=10,
        ),
        selection_strategy=selection_strategy.RouletteSelection(),
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
