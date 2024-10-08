from . import ActFuncs as ActFuncs
from typing import Any, Optional

class NetNode:
    def Eval(self, valVect: Any): ...
    inputNodes: Any = ...
    def SetInputs(self, inputNodes: Any) -> None: ...
    def GetInputs(self): ...
    weights: Any = ...
    def SetWeights(self, weights: Any) -> None: ...
    def GetWeights(self): ...
    nodeIndex: Any = ...
    nodeList: Any = ...
    actFunc: Any = ...
    def __init__(
        self,
        nodeIndex: Any,
        nodeList: Any,
        inputNodes: Optional[Any] = ...,
        weights: Optional[Any] = ...,
        actFunc: Any = ...,
        actFuncParms: Any = ...,
    ) -> None: ...
