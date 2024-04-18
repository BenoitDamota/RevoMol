from rdkit.ML.Data import cQuantize as cQuantize
from rdkit.ML.InfoTheory import entropy as entropy
from typing import Any

hascQuantize: int

def feq(v1: Any, v2: Any, tol: Any = ...): ...
def FindVarQuantBound(vals: Any, results: Any, nPossibleRes: Any): ...
def FindVarMultQuantBounds(
    vals: Any, nBounds: Any, results: Any, nPossibleRes: Any
): ...
