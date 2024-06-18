"""
Module to define the neighborhood strategy for a molecule representation.
"""

import random
from abc import ABC, abstractmethod

from evomol.action import (
    Action,
    pprint_action_space,
)
from evomol.representation import Molecule


class NeighborhoodStrategy(ABC):
    """
    Strategy to generate neighbors of a molecule representation.
    """

    depth: int

    @abstractmethod
    def mutate(self, molecule: Molecule) -> Molecule:
        """Takes a molecule and returns a new molecule with a mutation."""


class RandomNeighborhoodStrategy(NeighborhoodStrategy):
    """
    Random strategy to generate neighbors of a molecule representation.
    """

    def mutate(self, molecule: Molecule) -> Molecule:
        print(str(molecule))
        new_molecule: Molecule = Molecule(str(molecule))

        for depth in range(self.depth):
            print(f"Depth: {depth}")
            # list all possible actions
            possible_actions: dict[str, dict[str, list[Action]]] = (
                new_molecule.list_all_possible_actions()
            )

            if not possible_actions:
                raise ValueError("No possible actions")

            pprint_action_space(possible_actions)

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

            print(f"Applying action: {action}")
            new_molecule = action.apply()
            print(f"New molecule: {new_molecule}")

        return new_molecule
