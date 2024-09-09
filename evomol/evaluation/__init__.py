"""
This module contains the evaluation classes that are used to evaluate the
properties of the molecules.
"""

from .cl_score import CLScore as CLScore
from .cycle_score import CycleScore as CycleScore
from .cycle_score import NormalizedCycleScore as NormalizedCycleScore
from .diversity import Descriptor as Descriptor
from .diversity import Diversity as Diversity
from .diversity import atoms_list as atoms_list
from .diversity import checkmol as checkmol
from .diversity import compute_ifg as compute_ifg
from .diversity import ecfp4 as ecfp4
from .diversity import gen_scaffolds as gen_scaffolds
from .diversity import scaffolds as scaffolds
from .diversity import shg_1 as shg_1

# flake8: noqa
from .evaluation import Evaluation as Evaluation
from .evaluation import EvaluationError as EvaluationError
from .evaluation import Function as Function
from .evaluation import is_valid_molecule as is_valid_molecule
from .generic_cyclic_features import FilterUnknownGCF as FilterUnknownGCF
from .generic_cyclic_features import UnknownGCF as UnknownGCF
from .logp import LogP as LogP
from .logp import ZincNormalizedLogP as ZincNormalizedLogP
from .plogp import PLogP as PLogP
from .qed import QED as QED
from .rd_filters import RDFilters as RDFilters
from .sa_score import NormalizedSAScore as NormalizedSAScore
from .sa_score import SAScore as SAScore
from .sa_score import ZincNormalizedSAScore as ZincNormalizedSAScore
from .unknown_ecfp import FilterUnknownECFP as FilterUnknownECFP
from .unknown_ecfp import UnknownECFP as UnknownECFP

# removed to avoid guacamol dependency
# from .isomer import Isomer as Isomer
# from .rediscovery import Rediscovery as Rediscovery
