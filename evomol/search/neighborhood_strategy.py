"""
Module to define the neighborhood strategy for a molecule representation.
"""

import random
from abc import ABC, abstractmethod

from evomol.action import Action, ActionError
from evomol.evaluation import Evaluation, EvaluationError
from evomol.representation import Molecule
from evomol.search.parameters import search_parameters


class NoImproverError(RuntimeError):
    """Error raised when no neighbor is better than the current molecule."""


class NeighborhoodStrategy(ABC):
    """
    Strategy to generate neighbors of a molecule representation.
    """

    depth: int = 1
    nb_max_tries: int = 50

    @abstractmethod
    def mutate(
        self,
        molecule: Molecule,
        to_replace: Molecule,
        evaluations: list[Evaluation],
    ) -> list[Molecule]:
        """Takes a molecule and returns a new molecule with a mutation."""


class RandomNeighborhoodStrategy(NeighborhoodStrategy):
    """
    Random strategy to generate neighbors of a molecule representation.
    """

    def mutate(
        self,
        molecule: Molecule,
        to_replace: Molecule,
        evaluations: list[Evaluation],
    ) -> list[Molecule]:

        nb_attempts: int = 0
        while nb_attempts < self.nb_max_tries:
            if molecule.list_all_possible_actions() is None:
                break

            valid: bool = True

            # apply a random mutation
            try:
                new_molecule: Molecule = self.random_mutation(molecule)
            except ActionError as e:
                # print(e)
                valid = False
                raise e  # TODO ? should be a critical error ?

            for evaluation in evaluations:
                if not valid:
                    break
                try:
                    evaluation.evaluate(new_molecule)
                except EvaluationError:
                    # print(e)
                    valid = False

            # check if the new molecule is better than one of the molecules to replace
            if valid and (
                not to_replace
                or new_molecule.value(search_parameters.fitness_criteria)
                > to_replace.value(search_parameters.fitness_criteria)
            ):
                return new_molecule

            nb_attempts += 1

        return None

    def random_mutation(self, new_molecule: Molecule) -> Molecule:
        """Choose a random mutation and apply it to the molecule.

        Args:
            new_molecule (Molecule): Molecule to mutate.

        Returns:
            Molecule: New molecule after mutation.
        """

        for _ in range(self.depth):
            # print(f"Depth: {depth}")
            # list all possible actions
            possible_actions: dict[str, dict[str, list[Action]]] = (
                new_molecule.list_all_possible_actions()
            )

            # pprint_action_space(possible_actions)

            # select a random representation
            representation: str = random.choice(list(possible_actions.keys()))

            # select a random action space
            action_space: str = random.choice(
                list(possible_actions[representation].keys())
            )

            # select a random action
            action: Action = random.choice(
                possible_actions[representation][action_space]
            )

            # print(f"Applying action: {action}")
            new_molecule = action.apply()
            # print(f"New molecule: {new_molecule}")

        return new_molecule
