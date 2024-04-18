from rdkit import RDConfig as RDConfig
from typing import Any

def GetAtomicData(
    atomDict: Any,
    descriptorsDesired: Any,
    dBase: Any = ...,
    table: str = ...,
    where: str = ...,
    user: str = ...,
    password: str = ...,
    includeElCounts: int = ...,
) -> None: ...
def SplitComposition(compStr: Any): ...
def ConfigToNumElectrons(
    config: Any, ignoreFullD: int = ..., ignoreFullF: int = ...
): ...
