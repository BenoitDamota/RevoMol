from rdkit.VLib.Node import VLibNode as VLibNode
from typing import Any, Optional

class OutputNode(VLibNode):
    def __init__(
        self, dest: Optional[Any] = ..., strFunc: Optional[Any] = ..., **kwargs: Any
    ) -> None: ...
    def next(self): ...
