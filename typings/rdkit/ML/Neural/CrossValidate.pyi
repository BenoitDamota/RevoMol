from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.Neural import Network as Network, Trainers as Trainers
from typing import Any, Optional

def CrossValidate(
    net: Any, testExamples: Any, tolerance: Any, appendExamples: int = ...
): ...
def CrossValidationDriver(
    examples: Any,
    attrs: Any = ...,
    nPossibleVals: Any = ...,
    holdOutFrac: float = ...,
    silent: int = ...,
    tolerance: float = ...,
    calcTotalError: int = ...,
    hiddenSizes: Optional[Any] = ...,
    **kwargs: Any,
): ...
