"""
This module contains the classes for the molecular graph actions.
"""

from .add_atom import (
    ActionSpaceAddAtomMolGraph as ActionSpaceAddAtomMolGraph,
)
from .add_atom import (
    AddAtomMolGraph as AddAtomMolGraph,
)
from .change_bond import (
    ActionSpaceChangeBondMolGraph as ActionSpaceChangeBondMolGraph,
)
from .change_bond import (
    ChangeBondMolGraph as ChangeBondMolGraph,
)
from .cut_atom import (
    ActionSpaceCutAtomMolGraph as ActionSpaceCutAtomMolGraph,
)
from .cut_atom import (
    CutAtomMolGraph as CutAtomMolGraph,
)
from .insert_carbon import (
    ActionSpaceInsertCarbonMolGraph as ActionSpaceInsertCarbonMolGraph,
)
from .insert_carbon import (
    InsertCarbonMolGraph as InsertCarbonMolGraph,
)
from .move_functional_group import (
    ActionSpaceMoveFunctionalGroupMolGraph as ActionSpaceMoveFunctionalGroupMolGraph,
)
from .move_functional_group import (
    MoveFunctionalGroupMolGraph as MoveFunctionalGroupMolGraph,
)
from .remove_atom import (
    ActionSpaceRemoveAtomMolGraph as ActionSpaceRemoveAtomMolGraph,
)
from .remove_atom import (
    RemoveAtomMolGraph as RemoveAtomMolGraph,
)
from .remove_group import (
    ActionSpaceRemoveGroupMolGraph as ActionSpaceRemoveGroupMolGraph,
)
from .remove_group import (
    RemoveGroupMolGraph as RemoveGroupMolGraph,
)
from .substitute_atom import (
    ActionSpaceSubstituteAtomMolGraph as ActionSpaceSubstituteAtomMolGraph,
)
from .substitute_atom import (
    SubstituteAtomMolGraph as SubstituteAtomMolGraph,
)


# from evomol.representation.molecular_graph_actions.add_atom import *
# from evomol.representation.molecular_graph_actions.change_bond import *
# from .cut_atom import *
# from .insert_carbon import *
# from .move_functional_group import *
# from .remove_atom import *
# from .remove_group import *
# from .substitute_atom import *
