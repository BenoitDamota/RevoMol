from rdkit.ML.DecTree import DecTree as DecTree
from rdkit.ML.InfoTheory import entropy as entropy
from typing import Any, Optional

def CalcTotalEntropy(examples: Any, nPossibleVals: Any): ...
def GenVarTable(examples: Any, nPossibleVals: Any, vars: Any): ...
def ID3(
    examples: Any,
    target: Any,
    attrs: Any,
    nPossibleVals: Any,
    depth: int = ...,
    maxDepth: int = ...,
    **kwargs: Any,
): ...
def ID3Boot(
    examples: Any,
    attrs: Any,
    nPossibleVals: Any,
    initialVar: Optional[Any] = ...,
    depth: int = ...,
    maxDepth: int = ...,
    **kwargs: Any,
): ...
