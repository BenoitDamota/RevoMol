from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.KNN import DistFunctions as DistFunctions
from rdkit.ML.KNN.KNNClassificationModel import (
    KNNClassificationModel as KNNClassificationModel,
)
from rdkit.ML.KNN.KNNRegressionModel import KNNRegressionModel as KNNRegressionModel
from typing import Any

def makeClassificationModel(numNeigh: Any, attrs: Any, distFunc: Any): ...
def makeRegressionModel(numNeigh: Any, attrs: Any, distFunc: Any): ...
def CrossValidate(knnMod: Any, testExamples: Any, appendExamples: int = ...): ...
def CrossValidationDriver(
    examples: Any,
    attrs: Any,
    nPossibleValues: Any,
    numNeigh: Any,
    modelBuilder: Any = ...,
    distFunc: Any = ...,
    holdOutFrac: float = ...,
    silent: int = ...,
    calcTotalError: int = ...,
    **kwargs: Any,
): ...
