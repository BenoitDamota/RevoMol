from rdkit import RDRandom as RDRandom
from typing import Any

SeqTypes: Any

def SplitIndices(
    nPts: Any, frac: Any, silent: int = ..., legacy: int = ..., replacement: int = ...
): ...
def SplitDataSet(data: Any, frac: Any, silent: int = ...): ...
def SplitDbData(
    conn: Any,
    fracs: Any,
    table: str = ...,
    fields: str = ...,
    where: str = ...,
    join: str = ...,
    labelCol: str = ...,
    useActs: int = ...,
    nActs: int = ...,
    actCol: str = ...,
    actBounds: Any = ...,
    silent: int = ...,
): ...
