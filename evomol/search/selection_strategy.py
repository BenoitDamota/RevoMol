"""
This module contains the selection strategies for the evolutionary algorithm.

The selection strategy is used to select the order to mutate the molecules
in the population.
"""

import random
from abc import ABC, abstractmethod

from typing_extensions import override

from evomol.representation import Molecule
from evomol.search.parameters import search_parameters


class SelectionStrategy(ABC):
    """Abstract class for the selection strategy."""

    @abstractmethod
    def select(self, population: list[Molecule]) -> list[int]:
        """
        Select the molecules to be mutated.

        Args:
            population (list[Molecule]): The sorted population of molecules.

        Returns:
            list[int]: The order to select the molecules.
        """


class RandomSelection(SelectionStrategy):
    """Molecules are selected randomly."""

    @override
    def select(self, population: list[Molecule]) -> list[int]:
        """
        Select the molecules to be mutated.

        Args:
            population (list[Molecule]): The sorted population of molecules.

        Returns:
            list[int]: The order to select the molecules.
        """
        order = list(range(len(population)))
        random.shuffle(order)
        return order


class BestSelection(SelectionStrategy):
    """Best molecules are selected first."""

    @override
    def select(self, population: list[Molecule]) -> list[int]:
        """
        Select the molecules to be mutated.

        Args:
            population (list[Molecule]): The sorted population of molecules.

        Returns:
            list[int]: The order to select the molecules.
        """
        return list(range(len(population)))


class RouletteSelection(SelectionStrategy):
    """Best molecules have a higher chance to be selected first."""

    @override
    def select(self, population: list[Molecule]) -> list[int]:
        """
        Select the molecules to be mutated.

        Args:
            population (list[Molecule]): The sorted population of molecules.

        Returns:
            list[int]: The order to select the molecules.
        """
        # shuffle the order of the population but gives more chance to
        # the best to be at the beginning of the list

        total_fitness = sum(
            molecule.value(search_parameters.fitness_criteria)
            for molecule in population
        )

        selection_probabilities = [
            molecule.value(search_parameters.fitness_criteria) / total_fitness
            for molecule in population
        ]

        order_selection = random.choices(
            range(len(population)),
            weights=selection_probabilities,
            k=len(population),
        )
        return order_selection
