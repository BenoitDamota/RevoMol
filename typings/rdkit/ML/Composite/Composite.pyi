from rdkit.ML.Data import DataUtils as DataUtils
from typing import Any, Optional

class Composite:
    modelList: Any = ...
    errList: Any = ...
    countList: Any = ...
    modelVotes: Any = ...
    quantBounds: Any = ...
    nPossibleVals: Any = ...
    quantizationRequirements: Any = ...
    activityQuant: Any = ...
    def __init__(self) -> None: ...
    def SetModelFilterData(
        self, modelFilterFrac: float = ..., modelFilterVal: float = ...
    ) -> None: ...
    def SetDescriptorNames(self, names: Any) -> None: ...
    def GetDescriptorNames(self): ...
    def SetQuantBounds(self, qBounds: Any, nPossible: Optional[Any] = ...) -> None: ...
    def GetQuantBounds(self): ...
    def GetActivityQuantBounds(self): ...
    def SetActivityQuantBounds(self, bounds: Any) -> None: ...
    def QuantizeActivity(
        self, example: Any, activityQuant: Optional[Any] = ..., actCol: int = ...
    ): ...
    def QuantizeExample(self, example: Any, quantBounds: Optional[Any] = ...): ...
    def MakeHistogram(self): ...
    def CollectVotes(
        self,
        example: Any,
        quantExample: Any,
        appendExample: int = ...,
        onlyModels: Optional[Any] = ...,
    ): ...
    def ClassifyExample(
        self,
        example: Any,
        threshold: int = ...,
        appendExample: int = ...,
        onlyModels: Optional[Any] = ...,
    ): ...
    def GetVoteDetails(self): ...
    def GetInputOrder(self): ...
    def SetInputOrder(self, colNames: Any) -> None: ...
    def Grow(
        self,
        examples: Any,
        attrs: Any,
        nPossibleVals: Any,
        buildDriver: Any,
        pruner: Optional[Any] = ...,
        nTries: int = ...,
        pruneIt: int = ...,
        needsQuantization: int = ...,
        progressCallback: Optional[Any] = ...,
        **buildArgs: Any,
    ) -> None: ...
    def ClearModelExamples(self) -> None: ...
    def Pickle(self, fileName: str = ..., saveExamples: int = ...) -> None: ...
    def AddModel(
        self, model: Any, error: Any, needsQuantization: int = ...
    ) -> None: ...
    def AverageErrors(self): ...
    def SortModels(self, sortOnError: bool = ...) -> None: ...
    def GetModel(self, i: Any): ...
    def SetModel(self, i: Any, val: Any) -> None: ...
    def GetCount(self, i: Any): ...
    def SetCount(self, i: Any, val: Any) -> None: ...
    def GetError(self, i: Any): ...
    def SetError(self, i: Any, val: Any) -> None: ...
    def GetDataTuple(self, i: Any): ...
    def SetDataTuple(self, i: Any, tup: Any) -> None: ...
    def GetAllData(self): ...
    def __len__(self): ...
    def __getitem__(self, which: Any): ...
