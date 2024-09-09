"""
Module for the evolutionary algorithm.
"""

import time
from dataclasses import dataclass

from evomol.evaluation import Evaluation, is_valid_molecule
from evomol.logging import get_logger
from evomol.representation import MolecularGraph, Molecule
from evomol.search.neighborhood_strategy import NeighborhoodStrategy
from evomol.search.parameters import search_parameters
from evomol.search.selection_strategy import SelectionStrategy


@dataclass
class ParametersEvolutionaryAlgorithm:
    """Parameters for the evolutionary algorithm."""

    size: int = 1000
    # stop conditions
    max_generations: int = 1000
    max_time: int = 3600
    nb_replacements: int = 10


# pylint: disable=too-many-instance-attributes
class EvolutionaryAlgorithm:
    """
    Base class for evolutionary algorithms.
    """

    def __init__(
        self,
        population: list[Molecule],
        parameters: ParametersEvolutionaryAlgorithm,
        selection_strategy: SelectionStrategy,
        evaluations: list[Evaluation],
        neighborhood_strategy: NeighborhoodStrategy,
    ):
        self.population: list[Molecule] = population
        self.parameters: ParametersEvolutionaryAlgorithm = parameters
        self.selection_strategy: SelectionStrategy = selection_strategy
        self.evaluations: list[Evaluation] = evaluations
        self.neighborhood_strategy: NeighborhoodStrategy = neighborhood_strategy

        self.tabu_molecules: set[str] = set()

        if len(self.population) < self.parameters.size:
            self.population = self.population + [
                Molecule("") for _ in range(self.parameters.size - len(self.population))
            ]

        for molecule in self.population:
            self.tabu_molecules.add(
                molecule.get_representation(MolecularGraph).canonical_smiles
            )
            if not is_valid_molecule(molecule, evaluations):
                get_logger().error("Invalid molecule in the population : %s", molecule)
                raise ValueError("Invalid molecule in the population")

        # stop conditions
        self.nb_generations: int = 0
        self.start_time: float = 0
        self.stop_time: float = 0

    def run(self) -> None:
        """
        Run the evolutionary algorithm.
        """
        logger = get_logger()
        self.start_time = time.time()
        self.stop_time = self.start_time + self.parameters.max_time

        index_selected: int = 0

        logger.info("Start the evolutionary algorithm.")

        while self.continue_search():

            # sort the molecules for the selection
            # the order is used to select the molecules to mutate
            # that will replace the worst molecules
            order_selection = self.selection_strategy.select(self.population)

            # the first molecule of the order_selection list is the first
            # one to mutate
            current_selected = 0

            # for each of the worst molecules
            for index_to_replace in range(
                len(self.population) - 1,
                len(self.population) - 1 - self.parameters.nb_replacements,
                -1,
            ):

                found_improver: bool = False

                # while no improver is found and there are still molecules to
                # select
                while current_selected < len(order_selection) and not found_improver:

                    # select the molecule to mutate
                    index_selected = order_selection[current_selected]

                    # mutate the molecule
                    new_molecule: Molecule | None = self.neighborhood_strategy.mutate(
                        self.population[index_selected],
                        self.population[index_to_replace],
                        evaluations=self.evaluations,
                        tabu_molecules=self.tabu_molecules,
                    )

                    smiles_to_mutate = (
                        self.population[index_selected]
                        .get_representation(MolecularGraph)
                        .canonical_smiles
                    )
                    smiles_to_replace = (
                        self.population[index_to_replace]
                        .get_representation(MolecularGraph)
                        .canonical_smiles
                    )

                    # if a new molecule is found
                    if new_molecule is not None:
                        new_smiles = new_molecule.get_representation(
                            MolecularGraph
                        ).canonical_smiles
                        score = new_molecule.value(search_parameters.fitness_criteria)
                        logger.info(
                            "New molecule : %s (%1.4f) from %s to replace %s",
                            new_smiles,
                            score,
                            smiles_to_mutate,
                            smiles_to_replace,
                            extra={
                                "new_molecule": new_smiles,
                                "to_mutate": smiles_to_mutate,
                                "to_replace": smiles_to_replace,
                            },
                        )
                        # replace the molecule to replace
                        self.population[index_to_replace] = new_molecule
                        # add the new molecule to the tabu list
                        self.tabu_molecules.add(
                            new_molecule.get_representation(
                                MolecularGraph
                            ).canonical_smiles
                        )
                        found_improver = True
                    else:
                        logger.info(
                            "No improvement found from %s to replace %s",
                            smiles_to_mutate,
                            smiles_to_replace,
                            extra={
                                "to_mutate": smiles_to_mutate,
                                "to_replace": smiles_to_replace,
                            },
                        )

                    # go to the next molecule to select
                    current_selected += 1

            # sort the population by fitness
            self.population.sort(
                key=lambda x: x.value(search_parameters.fitness_criteria),
                reverse=True,
            )

            self.nb_generations += 1

            # statistics
            mean_value = sum(
                molecule.value(search_parameters.fitness_criteria)
                for molecule in self.population
            ) / len(self.population)

            max_value = self.population[0].value(search_parameters.fitness_criteria)
            min_value = self.population[-1].value(search_parameters.fitness_criteria)

            best_smiles = (
                self.population[0].get_representation(MolecularGraph).canonical_smiles
            )

            logger.info(
                "Generation %6d - Mean: %1.4f - Max: %1.4f - Min: %1.4f - %s",
                self.nb_generations,
                mean_value,
                max_value,
                min_value,
                best_smiles,
                extra={
                    "generation": self.nb_generations,
                    "mean": mean_value,
                    "max": max_value,
                    "min": min_value,
                    "best_smiles": best_smiles,
                },
            )

        elapsed_time = time.time() - self.start_time
        logger.info(
            "End of search - elapsed time: %s",
            elapsed_time,
            extra={"elapsed_time": elapsed_time},
        )
        bests = [
            (
                molecule.get_representation(MolecularGraph).canonical_smiles,
                molecule.value(search_parameters.fitness_criteria),
            )
            for molecule in self.population[:5]
        ]
        logger.info(
            "Bests: %s",
            bests,
            extra={"bests": {str(i): str(best) for i, best in enumerate(bests)}},
        )

    def continue_search(self) -> bool:
        """
        Check if the search should continue.
        """
        return (
            self.nb_generations < self.parameters.max_generations
            and self.stop_time > time.time()
        )
