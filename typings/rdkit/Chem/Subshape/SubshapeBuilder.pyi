from rdkit import Chem as Chem, Geometry as Geometry
from rdkit.Chem import AllChem as AllChem
from rdkit.Chem.Subshape import (
    BuilderUtils as BuilderUtils,
    SubshapeObjects as SubshapeObjects,
)
from typing import Any, Optional

class SubshapeCombineOperations:
    UNION: int = ...
    SUM: int = ...
    INTERSECT: int = ...

class SubshapeBuilder:
    gridDims: Any = ...
    gridSpacing: float = ...
    winRad: float = ...
    nbrCount: int = ...
    terminalPtRadScale: float = ...
    fraction: float = ...
    stepSize: float = ...
    featFactory: Any = ...
    def SampleSubshape(self, subshape1: Any, newSpacing: Any): ...
    def GenerateSubshapeShape(
        self, cmpd: Any, confId: int = ..., addSkeleton: bool = ..., **kwargs: Any
    ): ...
    def __call__(self, cmpd: Any, **kwargs: Any): ...
    def GenerateSubshapeSkeleton(
        self,
        shape: Any,
        conf: Optional[Any] = ...,
        terminalPtsOnly: bool = ...,
        skelFromConf: bool = ...,
    ) -> None: ...
    def CombineSubshapes(
        self, subshape1: Any, subshape2: Any, operation: Any = ...
    ): ...
