from typing import Any, Optional

class DescriptorCalculationError(Exception): ...
class ClassificationError(Exception): ...

class ModelPackage:
    def __init__(
        self,
        descCalc: Optional[Any] = ...,
        model: Optional[Any] = ...,
        dataSet: Optional[Any] = ...,
        notes: str = ...,
    ) -> None: ...
    def SetCalculator(self, calc: Any) -> None: ...
    def GetCalculator(self): ...
    def SetModel(self, model: Any) -> None: ...
    def GetModel(self): ...
    def SetDataset(self, data: Any) -> None: ...
    def GetDataset(self): ...
    def SetNotes(self, notes: Any) -> None: ...
    def GetNotes(self): ...
    def SetSupplementalData(self, suppD: Any) -> None: ...
    def GetSupplementalData(self): ...
    def AddSupplementalData(self, data: Any) -> None: ...
    def Classify(self, obj: Any, label: str = ..., threshold: int = ...): ...
    def Init(self) -> None: ...
