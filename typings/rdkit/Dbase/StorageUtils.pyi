from rdkit import RDConfig as RDConfig
from rdkit.Dbase import DbModule as DbModule
from typing import Any, Optional

def ValidateRDId(ID: Any): ...
def RDIdToInt(ID: Any, validate: int = ...): ...
def IndexToRDId(idx: Any, leadText: str = ...): ...
def GetNextId(conn: Any, table: Any, idColName: str = ...): ...
def GetNextRDId(conn: Any, table: Any, idColName: str = ..., leadText: str = ...): ...
def RegisterItem(
    conn: Any,
    table: Any,
    value: Any,
    columnName: Any,
    data: Optional[Any] = ...,
    id: str = ...,
    idColName: str = ...,
    leadText: str = ...,
): ...
def RegisterItems(
    conn: Any,
    table: Any,
    values: Any,
    columnName: Any,
    rows: Any,
    startId: str = ...,
    idColName: str = ...,
    leadText: str = ...,
): ...

__test__: Any
