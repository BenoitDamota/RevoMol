from rdkit.VLib.Node import VLibNode as VLibNode
from typing import Any, Optional

class SupplyNode(VLibNode):
    def __init__(self, contents: Optional[Any] = ..., **kwargs: Any) -> None: ...
    def reset(self) -> None: ...
    def next(self): ...
    def AddParent(self, parent: Any, notify: int = ...) -> None: ...