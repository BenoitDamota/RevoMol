"""
This module contains the classes for the molecular graph actions.
"""

# flake8: noqa
from .add_atom import (
    AddAtomMG as AddAtomMG,
)
from .add_group import (
    AddGroupMG as AddGroupMG,
)
from .change_bond import (
    ChangeBondMG as ChangeBondMG,
)
from .cut_atom import (
    CutAtomMG as CutAtomMG,
)
from .insert_carbon import (
    InsertCarbonMG as InsertCarbonMG,
)
from .move_group import (
    MoveGroupMG as MoveGroupMG,
)
from .remove_atom import (
    RemoveAtomMG as RemoveAtomMG,
)
from .remove_group import (
    RemoveGroupMG as RemoveGroupMG,
)
from .substitute_atom import (
    SubstituteAtomMG as SubstituteAtomMG,
)
