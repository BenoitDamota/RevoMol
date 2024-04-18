from rdkit.Chem import ChemicalFeatures as ChemicalFeatures
from typing import Any

class FeatMapPoint(ChemicalFeatures.FreeChemicalFeature):
    weight: float = ...
    featDirs: Any = ...
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    def initFromFeat(self, feat: Any) -> None: ...
    def GetDist2(self, other: Any): ...
    def GetDirMatch(self, other: Any, useBest: bool = ...): ...
