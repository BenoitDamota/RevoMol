"""
Module to define the neighborhood strategy for a molecule representation.
"""

import random
import sys
from abc import ABC, abstractmethod

from typing_extensions import override

import evomol.action.action as action_module
from evomol.action import Action, ActionError
from evomol.evaluation import Evaluation, is_valid_molecule
from evomol.representation import MolecularGraph, Molecule
from evomol.search.parameters import search_parameters


class NoImproverError(RuntimeError):
    """Error raised when no neighbor is better than the current molecule."""


class NeighborhoodStrategy(ABC):
    """
    Strategy to generate neighbors of a molecule representation.
    """

    def __init__(self, depth: int = 1, nb_max_tries: int = 50):
        self.depth = depth
        self.nb_max_tries = nb_max_tries

    @abstractmethod
    def mutate(
        self,
        molecule: Molecule,
        to_replace: Molecule,
        evaluations: list[Evaluation],
        tabu_molecules: set[str],
    ) -> Molecule | None:
        """Takes a molecule and returns a new molecule with a mutation."""


class RandomNeighborhoodStrategy(NeighborhoodStrategy):
    """
    Random strategy to generate neighbors of a molecule representation.
    """

    def __init__(
        self,
        depth: int = 1,
        nb_max_tries: int = 50,
        preselect_actions: bool = True,
    ):
        super().__init__(depth, nb_max_tries)
        # whether if one action space is selected before or if all actions
        # space are computed then randomly selected
        self.preselect_actions: bool = preselect_actions

    @override
    def mutate(
        self,
        molecule: Molecule,
        to_replace: Molecule,
        evaluations: list[Evaluation],
        tabu_molecules: set[str],
    ) -> Molecule | None:

        nb_attempts: int = 0
        while nb_attempts < self.nb_max_tries:

            # apply a random mutation
            try:
                new_molecule: Molecule = self.random_mutation(molecule)
            except ActionError as e:
                # Is this a critical error? This would mean that the action
                # is invalid and the code is not working properly.
                # as only update_representation is supposed to crash
                print(e, file=sys.stderr)
                # count as a failed mutation
                nb_attempts += 1
                continue

            new_smiles: str = new_molecule.get_representation(
                MolecularGraph
            ).canonical_smiles
            # check if the new molecule is in the tabu list
            if new_smiles in tabu_molecules:
                nb_attempts += 1
                continue

            evaluation_ok = is_valid_molecule(new_molecule, evaluations)

            # check if the new molecule is better than one of the molecules to replace
            if evaluation_ok and (
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
            action: Action
            if self.preselect_actions:
                # select the action space before generating the actions

                representation_dict: dict[str, list[type[action_module.Action]]] = (
                    new_molecule.list_available_actions_space()
                )

                # choose a random representation
                representation_choice: str = random.choice(
                    list(representation_dict.keys())
                )

                # choose a random action
                action_choice: type[action_module.Action] = random.choice(
                    representation_dict[representation_choice]
                )

                # list actions
                actions: list[action_module.Action] = action_choice.list_actions(
                    new_molecule
                )

                if not actions:
                    continue

                # choose a random action
                action = random.choice(actions)

            else:
                # generate all possible actions then choose one randomly
                possible_actions: dict[str, dict[str, list[action_module.Action]]] = (
                    new_molecule.compute_possible_actions()
                )

                # if there are no possible actions, return the molecule
                if not possible_actions:
                    return new_molecule

                # choose a random representation
                representation: str = random.choice(list(possible_actions.keys()))

                # choose a random action space
                action_space: str = random.choice(
                    list(possible_actions[representation].keys())
                )

                # choose a random action
                action = random.choice(possible_actions[representation][action_space])

            # apply the action
            new_molecule = action.apply()
        return new_molecule
