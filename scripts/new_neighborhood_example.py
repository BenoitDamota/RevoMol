"""
This script launches an evolutionary algorithm to maximize the QED of a molecule.
"""

import os
import random
import sys

from typing_extensions import override

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error, R0801
# noqa: E402

from evomol import default_parameters as dp
from evomol import evaluation as evaluator
from evomol.logging import init_logger
from evomol.representation import MolecularGraph, Molecule
from evomol.search import selection_strategy
from evomol.search.enumeration import find_neighbors
from evomol.search.evolutionary_algorithm import (
    EvolutionaryAlgorithm,
    ParametersEvolutionaryAlgorithm,
)
from evomol.search.neighborhood_strategy import NeighborhoodStrategy
from evomol.search.parameters import search_parameters


class EnumerationNeighborhoodStrategy(NeighborhoodStrategy):
    """
    Strategy that consists in enumerating all the neighbors of a molecule, filter
    them and save them.
    """

    def __init__(self, depth: int, nb_max_tries: int):
        super().__init__(depth, nb_max_tries)
        self.found_molecules: dict[str, list[str]] = {}

    @override
    def mutate(
        self,
        molecule: Molecule,
        to_replace: Molecule,
        evaluations: list[evaluator.Evaluation],
        tabu_molecules: set[str],
    ) -> Molecule | None:
        smiles = molecule.get_representation(MolecularGraph).canonical_smiles

        # if the molecule has not been explored yet, we need to find its neighbors
        if smiles not in self.found_molecules:
            self.found_molecules[smiles] = []
            neighbors: set[str] = find_neighbors(molecule, self.depth)
            # check if the neighbors are valid
            for neighbor in neighbors:
                is_valid = True
                neighbor_ = Molecule(neighbor)
                if neighbor in tabu_molecules:
                    continue
                for eval_ in evaluations:
                    # if the eval is unknown ecfp or gcf, we need to evaluate the
                    # molecule first
                    if "UnknownECFP" in eval_.name or "UnknownGCF" in eval_.name:
                        try:
                            eval_.evaluate(neighbor_)
                        except evaluator.EvaluationError:
                            is_valid = False
                            break
                if not is_valid:
                    continue
                # save the neighbor
                self.found_molecules[smiles].append(neighbor)

        # if all the neighbors have been explored, return None
        if not self.found_molecules[smiles]:
            return None

        # try to find a neighbor with a good enough fitness
        nb_attempts = 0
        while nb_attempts < self.nb_max_tries and self.found_molecules[smiles][0]:
            nb_attempts += 1
            neighbor = random.choice(self.found_molecules[smiles])
            self.found_molecules[smiles].remove(neighbor)
            if neighbor in tabu_molecules:
                continue
            neighbor_ = Molecule(neighbor)
            evaluation_ok = evaluator.is_valid_molecule(neighbor_, evaluations)

            # check if the new molecule is better than one of the molecules to replace
            if evaluation_ok and (
                not to_replace
                or neighbor_.value(search_parameters.fitness_criteria)
                > to_replace.value(search_parameters.fitness_criteria)
            ):
                return neighbor_

        # if no neighbor has been found, return None
        return None


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
        config_file="logging_configs/stderr-json-file.json",
        output_file="logs/evomol_default.log",
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
            size=20,
            max_generations=1000,
            max_time=3600,
            nb_replacements=3,
        ),
        selection_strategy=selection_strategy.RouletteSelection(),
        evaluations=evaluations,
        neighborhood_strategy=EnumerationNeighborhoodStrategy(
            depth=2,
            nb_max_tries=50,
        ),
    )

    # launch the search
    evo_mol.run()


if __name__ == "__main__":
    main()
