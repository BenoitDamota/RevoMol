"""
This module contains the actions used on the molecule sorted by the type of
representation used.
Each action is a class that inherits from the Action class and implements the
necessary methods to perform the action and list possible actions on a molecule.
"""

# flake8: noqa
from .action import Action as Action
from .action import ActionError as ActionError
from .action import pprint_action_space as pprint_action_space
from .action import pprint_action_space_and_apply as pprint_action_space_and_apply
