from rdkit import Chem as Chem
from rdkit.VLib.Output import OutputNode as BaseOutputNode
from typing import Any, Optional

class OutputNode(BaseOutputNode):
    def __init__(
        self,
        dest: Optional[Any] = ...,
        delim: str = ...,
        idField: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def reset(self) -> None: ...
    def smilesOut(self, mol: Any): ...
