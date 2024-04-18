from rdkit import RDConfig as RDConfig
from rdkit.Dbase import DbModule as DbModule
from typing import Any, Optional

sqlTextTypes: Any
sqlIntTypes: Any
sqlFloatTypes: Any
sqlBinTypes: Any

def GetDbNames(
    user: str = ...,
    password: str = ...,
    dirName: str = ...,
    dBase: str = ...,
    cn: Optional[Any] = ...,
): ...
def GetTableNames(
    dBase: Any,
    user: str = ...,
    password: str = ...,
    includeViews: int = ...,
    cn: Optional[Any] = ...,
): ...
def GetColumnInfoFromCursor(cursor: Any): ...
def GetColumnNamesAndTypes(
    dBase: Any,
    table: Any,
    user: str = ...,
    password: str = ...,
    join: str = ...,
    what: str = ...,
    cn: Optional[Any] = ...,
): ...
def GetColumnNames(
    dBase: Any,
    table: Any,
    user: str = ...,
    password: str = ...,
    join: str = ...,
    what: str = ...,
    cn: Optional[Any] = ...,
): ...
