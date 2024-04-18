from rdkit import DataStructs as DataStructs, RDConfig as RDConfig
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from rdkit.ML import CompositeRun as CompositeRun
from rdkit.ML.Data import DataUtils as DataUtils, SplitData as SplitData, Stats as Stats
from typing import Any

def cmp(t1: Any, t2: Any): ...
def message(msg: Any, noRet: int = ..., dest: Any = ...) -> None: ...
def error(msg: Any, dest: Any = ...) -> None: ...
def ScreenModel(
    mdl: Any,
    descs: Any,
    data: Any,
    picking: Any = ...,
    indices: Any = ...,
    errorEstimate: int = ...,
): ...
def AccumulateCounts(predictions: Any, thresh: int = ..., sortIt: int = ...): ...
def MakePlot(
    details: Any,
    final: Any,
    counts: Any,
    pickVects: Any,
    nModels: Any,
    nTrueActs: int = ...,
) -> None: ...
def Usage() -> None: ...
