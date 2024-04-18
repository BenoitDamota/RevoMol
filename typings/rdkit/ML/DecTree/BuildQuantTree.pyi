from rdkit.ML.Data import Quantize as Quantize
from rdkit.ML.DecTree import ID3 as ID3, QuantTree as QuantTree
from rdkit.ML.InfoTheory import entropy as entropy
from typing import Any, Optional

def FindBest(
    resCodes: Any,
    examples: Any,
    nBoundsPerVar: Any,
    nPossibleRes: Any,
    nPossibleVals: Any,
    attrs: Any,
    exIndices: Optional[Any] = ...,
    **kwargs: Any,
): ...
def BuildQuantTree(
    examples: Any,
    target: Any,
    attrs: Any,
    nPossibleVals: Any,
    nBoundsPerVar: Any,
    depth: int = ...,
    maxDepth: int = ...,
    exIndices: Optional[Any] = ...,
    **kwargs: Any,
): ...
def QuantTreeBoot(
    examples: Any,
    attrs: Any,
    nPossibleVals: Any,
    nBoundsPerVar: Any,
    initialVar: Optional[Any] = ...,
    maxDepth: int = ...,
    **kwargs: Any,
): ...
def TestTree() -> None: ...
def TestQuantTree() -> None: ...
def TestQuantTree2() -> None: ...
