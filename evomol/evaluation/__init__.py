"""
This module contains the evaluation classes that are used to evaluate the
properties of the molecules.
"""

# flake8: noqa
from .evaluation import (
    EvaluationError as EvaluationError,
    Evaluation as Evaluation,
    Function as Function,
    is_valid_molecule as is_valid_molecule,
)

from .cl_score import CLScore as CLScore
from .cycle_score import (
    CycleScore as CycleScore,
    NormalizedCycleScore as NormalizedCycleScore,
)


from .generic_cyclic_features import (
    UnknownGCF as UnknownGCF,
    FilterUnknownGCF as FilterUnknownGCF,
)

from .logp import LogP as LogP, ZincNormalizedLogP as ZincNormalizedLogP
from .qed import QED as QED
from .rd_filters import RDFilters as RDFilters
from .sa_score import (
    SAScore as SAScore,
    NormalizedSAScore as NormalizedSAScore,
    ZincNormalizedSAScore as ZincNormalizedSAScore,
)
from .unknown_ecfp import (
    UnknownECFP as UnknownECFP,
    FilterUnknownECFP as FilterUnknownECFP,
)

from .diversity import (
    Descriptor as Descriptor,
    Diversity as Diversity,
    scaffolds as scaffolds,
    gen_scaffolds as gen_scaffolds,
    compute_ifg as compute_ifg,
    atoms_list as atoms_list,
    shg_1 as shg_1,
    checkmol as checkmol,
    ecfp4 as ecfp4,
)

# removed to avoid guacamol dependency
# from .isomer import Isomer as Isomer
# from .rediscovery import Rediscovery as Rediscovery
