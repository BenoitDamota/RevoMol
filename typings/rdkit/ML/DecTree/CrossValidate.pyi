from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.DecTree import ID3 as ID3, randomtest as randomtest
from typing import Any, Optional

def ChooseOptimalRoot(
    examples: Any,
    trainExamples: Any,
    testExamples: Any,
    attrs: Any,
    nPossibleVals: Any,
    treeBuilder: Any,
    nQuantBounds: Any = ...,
    **kwargs: Any,
): ...
def CrossValidate(tree: Any, testExamples: Any, appendExamples: int = ...): ...
def CrossValidationDriver(
    examples: Any,
    attrs: Any,
    nPossibleVals: Any,
    holdOutFrac: float = ...,
    silent: int = ...,
    calcTotalError: int = ...,
    treeBuilder: Any = ...,
    lessGreedy: int = ...,
    startAt: Optional[Any] = ...,
    nQuantBounds: Any = ...,
    maxDepth: int = ...,
    **kwargs: Any,
): ...
def TestRun() -> None: ...
