from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.FeatureSelect import CMIM as CMIM
from rdkit.ML.NaiveBayes.ClassificationModel import (
    NaiveBayesClassifier as NaiveBayesClassifier,
)
from typing import Any, Optional

def makeNBClassificationModel(
    trainExamples: Any,
    attrs: Any,
    nPossibleValues: Any,
    nQuantBounds: Any,
    mEstimateVal: Any = ...,
    useSigs: bool = ...,
    ensemble: Optional[Any] = ...,
    useCMIM: int = ...,
    **kwargs: Any,
): ...
def CrossValidate(NBmodel: Any, testExamples: Any, appendExamples: int = ...): ...
def CrossValidationDriver(
    examples: Any,
    attrs: Any,
    nPossibleValues: Any,
    nQuantBounds: Any,
    mEstimateVal: float = ...,
    holdOutFrac: float = ...,
    modelBuilder: Any = ...,
    silent: int = ...,
    calcTotalError: int = ...,
    **kwargs: Any,
): ...
