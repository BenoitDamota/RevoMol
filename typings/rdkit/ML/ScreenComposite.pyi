from rdkit import DataStructs as DataStructs
from rdkit.Dbase import DbModule as DbModule
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from rdkit.ML import CompositeRun as CompositeRun
from rdkit.ML.Data import DataUtils as DataUtils, SplitData as SplitData
from typing import Any, Optional

hasPil: int

def message(msg: Any, noRet: int = ...) -> None: ...
def error(msg: Any) -> None: ...
def CalcEnrichment(mat: Any, tgt: int = ...): ...
def CollectResults(
    indices: Any,
    dataSet: Any,
    composite: Any,
    callback: Optional[Any] = ...,
    appendExamples: int = ...,
    errorEstimate: int = ...,
): ...
def DetailedScreen(
    indices: Any,
    data: Any,
    composite: Any,
    threshold: int = ...,
    screenResults: Optional[Any] = ...,
    goodVotes: Optional[Any] = ...,
    badVotes: Optional[Any] = ...,
    noVotes: Optional[Any] = ...,
    callback: Optional[Any] = ...,
    appendExamples: int = ...,
    errorEstimate: int = ...,
) -> None: ...
def ShowVoteResults(
    indices: Any,
    data: Any,
    composite: Any,
    nResultCodes: Any,
    threshold: Any,
    verbose: int = ...,
    screenResults: Optional[Any] = ...,
    callback: Optional[Any] = ...,
    appendExamples: int = ...,
    goodVotes: Optional[Any] = ...,
    badVotes: Optional[Any] = ...,
    noVotes: Optional[Any] = ...,
    errorEstimate: int = ...,
): ...
def ScreenIt(
    composite: Any,
    indices: Any,
    data: Any,
    partialVote: int = ...,
    voteTol: float = ...,
    verbose: int = ...,
    screenResults: Optional[Any] = ...,
    goodVotes: Optional[Any] = ...,
    badVotes: Optional[Any] = ...,
    noVotes: Optional[Any] = ...,
): ...
def PrepareDataFromDetails(model: Any, details: Any, data: Any, verbose: int = ...): ...
def ScreenFromDetails(
    models: Any,
    details: Any,
    callback: Optional[Any] = ...,
    setup: Optional[Any] = ...,
    appendExamples: int = ...,
    goodVotes: Optional[Any] = ...,
    badVotes: Optional[Any] = ...,
    noVotes: Optional[Any] = ...,
    data: Optional[Any] = ...,
    enrichments: Optional[Any] = ...,
): ...
def GetScreenImage(nGood: Any, nBad: Any, nRej: Any, size: Optional[Any] = ...): ...
def ScreenToHtml(
    nGood: Any,
    nBad: Any,
    nRej: Any,
    avgGood: Any,
    avgBad: Any,
    avgSkip: Any,
    voteTable: Any,
    imgDir: str = ...,
    fullPage: int = ...,
    skipImg: int = ...,
    includeDefs: int = ...,
): ...
def MakePredPlot(
    details: Any,
    indices: Any,
    data: Any,
    goodVotes: Any,
    badVotes: Any,
    nRes: Any,
    idCol: int = ...,
    verbose: int = ...,
) -> None: ...
def Go(details: Any) -> None: ...
def SetDefaults(details: Optional[Any] = ...): ...
def Usage() -> None: ...
def ShowVersion(includeArgs: int = ...) -> None: ...
def ParseArgs(details: Any): ...
