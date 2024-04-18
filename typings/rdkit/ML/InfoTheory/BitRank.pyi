from rdkit.ML.InfoTheory import entropy as entropy
from typing import Any

def FormCounts(
    bitVects: Any,
    actVals: Any,
    whichBit: Any,
    nPossibleActs: Any,
    nPossibleBitVals: int = ...,
): ...
def CalcInfoGains(
    bitVects: Any, actVals: Any, nPossibleActs: Any, nPossibleBitVals: int = ...
): ...
def RankBits(
    bitVects: Any, actVals: Any, nPossibleBitVals: int = ..., metricFunc: Any = ...
): ...
def AnalyzeSparseVects(bitVects: Any, actVals: Any): ...
def SparseRankBits(bitVects: Any, actVals: Any, metricFunc: Any = ...): ...