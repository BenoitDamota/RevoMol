from rdkit import Chem as Chem
from rdkit.VLib.Supply import SupplyNode as SupplyNode
from typing import Any

class SmilesSupplyNode(SupplyNode):
    def __init__(
        self,
        fileName: Any,
        delim: str = ...,
        nameColumn: int = ...,
        smilesColumn: int = ...,
        titleLine: int = ...,
        **kwargs: Any,
    ) -> None: ...
    def reset(self) -> None: ...
    def next(self): ...
