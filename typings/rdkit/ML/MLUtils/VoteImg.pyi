from typing import Any

def CollectVotes(composite: Any, data: Any, badOnly: Any): ...
def BuildVoteImage(
    nModels: Any,
    data: Any,
    values: Any,
    trueValues: Any = ...,
    sortTrueVals: int = ...,
    xScale: int = ...,
    yScale: int = ...,
    addLine: int = ...,
): ...
def VoteAndBuildImage(
    composite: Any,
    data: Any,
    badOnly: int = ...,
    sortTrueVals: int = ...,
    xScale: int = ...,
    yScale: int = ...,
    addLine: int = ...,
): ...
def Usage() -> None: ...
