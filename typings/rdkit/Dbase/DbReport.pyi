from rdkit import Chem as Chem
from rdkit.Chem import rdDepictor as rdDepictor
from rdkit.Chem.Draw import DrawUtils as DrawUtils
from rdkit.Chem.Draw.MolDrawing import (
    DrawingOptions as DrawingOptions,
    MolDrawing as MolDrawing,
)
from rdkit.Dbase import DbInfo as DbInfo
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from rdkit.Reports.PDFReport import PDFReport as PDFReport, ReportUtils as ReportUtils
from rdkit.utils import cactvs as cactvs
from typing import Any, Optional

GetReportlabTable: Any
QuickReport: Any
hasCDX: int

class CDXImageTransformer:
    smiCol: Any = ...
    tempHandler: Any = ...
    width: Any = ...
    verbose: Any = ...
    def __init__(
        self,
        smiCol: Any,
        width: int = ...,
        verbose: int = ...,
        tempHandler: Optional[Any] = ...,
    ) -> None: ...
    def __call__(self, arg: Any): ...

class CactvsImageTransformer:
    smiCol: Any = ...
    tempHandler: Any = ...
    width: Any = ...
    verbose: Any = ...
    def __init__(
        self,
        smiCol: Any,
        width: float = ...,
        verbose: int = ...,
        tempHandler: Optional[Any] = ...,
    ) -> None: ...
    def __call__(self, arg: Any): ...

class ReportLabImageTransformer:
    smiCol: Any = ...
    width: Any = ...
    verbose: Any = ...
    def __init__(
        self,
        smiCol: Any,
        width: float = ...,
        verbose: int = ...,
        tempHandler: Optional[Any] = ...,
    ) -> None: ...
    def __call__(self, arg: Any): ...

class RDImageTransformer:
    smiCol: Any = ...
    tempHandler: Any = ...
    width: Any = ...
    verbose: Any = ...
    def __init__(
        self,
        smiCol: Any,
        width: float = ...,
        verbose: int = ...,
        tempHandler: Optional[Any] = ...,
    ) -> None: ...
    def __call__(self, arg: Any): ...
