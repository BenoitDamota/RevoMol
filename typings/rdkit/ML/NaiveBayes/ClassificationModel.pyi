from rdkit.ML.Data import Quantize as Quantize
from typing import Any

class NaiveBayesClassifier:
    mprob: Any = ...
    def __init__(
        self,
        attrs: Any,
        nPossibleVals: Any,
        nQuantBounds: Any,
        mEstimateVal: Any = ...,
        useSigs: bool = ...,
    ) -> None: ...
    def GetName(self): ...
    def SetName(self, name: Any) -> None: ...
    def NameModel(self, varNames: Any) -> None: ...
    def GetExamples(self): ...
    def SetExamples(self, examples: Any) -> None: ...
    def GetTrainingExamples(self): ...
    def SetTrainingExamples(self, examples: Any) -> None: ...
    def GetTestExamples(self): ...
    def SetTestExamples(self, examples: Any) -> None: ...
    def SetBadExamples(self, examples: Any) -> None: ...
    def GetBadExamples(self): ...
    def trainModel(self) -> None: ...
    def ClassifyExamples(self, examples: Any, appendExamples: int = ...): ...
    def GetClassificationDetails(self): ...
    def ClassifyExample(self, example: Any, appendExamples: int = ...): ...
