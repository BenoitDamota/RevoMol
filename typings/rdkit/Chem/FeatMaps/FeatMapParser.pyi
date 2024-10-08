from rdkit import Geometry as Geometry
from rdkit.Chem.FeatMaps import FeatMapPoint as FeatMapPoint, FeatMaps as FeatMaps
from typing import Any, Optional

class FeatMapParseError(ValueError): ...

class FeatMapParser:
    data: Any = ...
    def __init__(
        self, file: Optional[Any] = ..., data: Optional[Any] = ...
    ) -> None: ...
    def SetData(self, data: Any) -> None: ...
    def Parse(self, featMap: Optional[Any] = ...): ...
    def ParseParamBlock(self): ...
    def ParseFeatPointBlock(self): ...
