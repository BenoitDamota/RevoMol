from rdkit.Dbase import DbInfo as DbInfo, DbModule as DbModule
from rdkit.Dbase.DbResultSet import (
    DbResultSet as DbResultSet,
    RandomAccessDbResultSet as RandomAccessDbResultSet,
)
from typing import Any, Optional

def GetColumns(
    dBase: Any,
    table: Any,
    fieldString: Any,
    user: str = ...,
    password: str = ...,
    join: str = ...,
    cn: Optional[Any] = ...,
): ...
def GetData(
    dBase: Any,
    table: Any,
    fieldString: str = ...,
    whereString: str = ...,
    user: str = ...,
    password: str = ...,
    removeDups: int = ...,
    join: str = ...,
    forceList: int = ...,
    transform: Optional[Any] = ...,
    randomAccess: int = ...,
    extras: Optional[Any] = ...,
    cn: Optional[Any] = ...,
): ...
def DatabaseToText(
    dBase: Any,
    table: Any,
    fields: str = ...,
    join: str = ...,
    where: str = ...,
    user: str = ...,
    password: str = ...,
    delim: str = ...,
    cn: Optional[Any] = ...,
): ...
def TypeFinder(data: Any, nRows: Any, nCols: Any, nullMarker: Optional[Any] = ...): ...
def GetTypeStrings(colHeadings: Any, colTypes: Any, keyCol: Optional[Any] = ...): ...
def TextFileToDatabase(
    dBase: Any,
    table: Any,
    inF: Any,
    delim: str = ...,
    user: str = ...,
    password: str = ...,
    maxColLabelLen: int = ...,
    keyCol: Optional[Any] = ...,
    nullMarker: Optional[Any] = ...,
) -> None: ...
def DatabaseToDatabase(
    fromDb: Any,
    fromTbl: Any,
    toDb: Any,
    toTbl: Any,
    fields: str = ...,
    join: str = ...,
    where: str = ...,
    user: str = ...,
    password: str = ...,
    keyCol: Optional[Any] = ...,
    nullMarker: str = ...,
) -> None: ...
