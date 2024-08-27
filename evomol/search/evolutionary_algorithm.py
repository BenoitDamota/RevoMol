"""
Module for the evolutionary algorithm.
"""

import random
import time
from dataclasses import dataclass

from evomol.evaluation import Evaluation
from evomol.representation import MolecularGraph, Molecule
from evomol.search.neighborhood_strategy import NeighborhoodStrategy
from evomol.search.parameters import search_parameters


@dataclass
class ParametersEvolutionaryAlgorithm:
    """Parameters for the evolutionary algorithm."""

    size: int = 1000
    # stop conditions
    max_generations: int = 1000
    max_time: int = 3600
    # selection
    nb_replacements: int = 10
    selection: str = "random"  # random, best, roulette


# pylint: disable=too-many-instance-attributes
class EvolutionaryAlgorithm:
    """
    Base class for evolutionary algorithms.
    """

    def __init__(
        self,
        population: list[Molecule],
        parameters: ParametersEvolutionaryAlgorithm,
        evaluations: list[Evaluation],
        neighborhood_strategy: NeighborhoodStrategy,
    ):
        self.population: list[Molecule] = population
        self.parameters: ParametersEvolutionaryAlgorithm = parameters
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
                for evaluation in evaluations:
                    evaluation.evaluate(molecule)

        # stop conditions
        self.nb_generations: int = 0
        self.start_time: float = 0
        self.stop_time: float = 0

    def run(self) -> None:
        """
        Run the evolutionary algorithm.
        """
        self.start_time = time.time()
        self.stop_time = self.start_time + self.parameters.max_time

        index_selected: int = 0

        while self.continue_search():

            order_selection = list(range(len(self.population)))
            if self.parameters.selection == "best":
                pass  # already sorted
            elif self.parameters.selection == "random":
                random.shuffle(order_selection)
            elif self.parameters.selection == "roulette":
                # shuffle the order of the population but gives more chance to
                # the best to be at the beginning of the list

                total_fitness = sum(
                    molecule.value(search_parameters.fitness_criteria)
                    for molecule in self.population
                )

                selection_probabilities = [
                    molecule.value(search_parameters.fitness_criteria) / total_fitness
                    for molecule in self.population
                ]

                order_selection = random.choices(
                    range(len(self.population)),
                    weights=selection_probabilities,
                    k=len(self.population),
                )

            # index_selected = random.randint(0, min(len(self.population), 5) - 1)
            # nb_possible_actions =
            # self.population[index_selected].nb_possible_actions()
            # while nb_possible_actions == 0:
            #     index_selected = index_selected + 1 % len(self.population)
            #     nb_possible_actions = self.population[
            #         index_selected
            #     ].nb_possible_actions()

            # to_replace = (
            #     self.population[self.parameters.nb_replacements :]
            #     if len(self.population) > self.parameters.size
            #     else []
            # )

            current_selected = 0

            for index_to_replace in range(
                len(self.population) - 1,
                len(self.population) - 1 - self.parameters.nb_replacements,
                -1,
            ):
                molecule_to_replace = self.population[index_to_replace]

                found_improver: bool = False

                while current_selected < len(order_selection) and not found_improver:
                    index_selected = order_selection[current_selected]

                    new_molecule: Molecule | None = self.neighborhood_strategy.mutate(
                        self.population[index_selected],
                        molecule_to_replace,
                        evaluations=self.evaluations,
                        tabu_molecules=self.tabu_molecules,
                    )

                    if new_molecule is not None:
                        self.population[index_to_replace] = new_molecule
                        self.tabu_molecules.add(
                            new_molecule.get_representation(
                                MolecularGraph
                            ).canonical_smiles
                        )
                        found_improver = True

                    current_selected += 1

            self.population.sort(
                key=lambda x: x.value(search_parameters.fitness_criteria),
                reverse=True,
            )

            self.nb_generations += 1
            if self.nb_generations % 10 == 0:
                _ = self.nb_generations

            mean_value = sum(
                molecule.value(search_parameters.fitness_criteria)
                for molecule in self.population
            ) / len(self.population)

            max_value = self.population[0].value(search_parameters.fitness_criteria)
            min_value = self.population[-1].value(search_parameters.fitness_criteria)

            best_smiles = (
                self.population[0].get_representation(MolecularGraph).canonical_smiles
            )

            print(
                f"-- Generation {self.nb_generations:6} - "
                # f"Selected {index_selected:2} "
                # f"({nb_possible_actions:3}-{selected.nb_possible_actions():3}) - "
                # f"nb new {nb_new_molecules:2} - "
                # f"{selected.get_representation(MolecularGraph).smiles} \n"
                f"Mean: {mean_value:1.4f} - "
                f"Max: {max_value:1.4f} - "
                f"Min: {min_value:1.4f} - "
                f" {best_smiles}"
            )
        print("End of search")
        bests = [
            (
                molecule.get_representation(MolecularGraph).canonical_smiles,
                molecule.value(search_parameters.fitness_criteria),
            )
            for molecule in self.population[:5]
        ]
        print(f"Bests:{bests}")

    def continue_search(self) -> bool:
        """
        Check if the search should continue.
        """
        return (
            self.nb_generations < self.parameters.max_generations
            and self.stop_time > time.time()
        )
