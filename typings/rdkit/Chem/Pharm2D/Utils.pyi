from typing import Any, Optional

nPointDistDict: Any
nDistPointDict: Any

def GetTriangles(nPts: Any): ...
def BinsTriangleInequality(d1: Any, d2: Any, d3: Any): ...
def ScaffoldPasses(combo: Any, bins: Optional[Any] = ...): ...
def NumCombinations(nItems: Any, nSlots: Any): ...
def CountUpTo(
    nItems: Any, nSlots: Any, vs: Any, idx: int = ..., startAt: int = ...
): ...
def GetIndexCombinations(
    nItems: Any, nSlots: Any, slot: int = ..., lastItemVal: int = ...
): ...
def GetAllCombinations(choices: Any, noDups: int = ..., which: int = ...): ...
def GetUniqueCombinations(choices: Any, classes: Any, which: int = ...): ...
def GetUniqueCombinations_new(choices: Any, classes: Any, which: int = ...): ...
def UniquifyCombinations(combos: Any): ...
def GetPossibleScaffolds(nPts: Any, bins: Any, useTriangleInequality: bool = ...): ...
def OrderTriangle(featIndices: Any, dists: Any): ...
