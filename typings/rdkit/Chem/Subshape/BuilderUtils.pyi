from rdkit import Geometry as Geometry
from rdkit.Chem.Subshape import SubshapeObjects as SubshapeObjects
from typing import Any

def ComputeGridIndices(shapeGrid: Any, winRad: Any): ...
def ComputeShapeGridCentroid(pt: Any, shapeGrid: Any, winRad: Any): ...
def FindTerminalPtsFromShape(
    shape: Any, winRad: Any, fraction: Any, maxGridVal: int = ...
): ...
def FindTerminalPtsFromConformer(conf: Any, winRad: Any, nbrCount: Any): ...
def FindGridPointBetweenPoints(pt1: Any, pt2: Any, shapeGrid: Any, winRad: Any): ...
def ClusterTerminalPts(pts: Any, winRad: Any, scale: Any): ...
def GetMoreTerminalPoints(
    shape: Any, pts: Any, winRad: Any, maxGridVal: Any, targetNumber: int = ...
) -> None: ...
def FindFarthestGridPoint(shape: Any, loc: Any, winRad: Any, maxGridVal: Any): ...
def ExpandTerminalPts(
    shape: Any, pts: Any, winRad: Any, maxGridVal: float = ..., targetNumPts: int = ...
) -> None: ...
def AppendSkeletonPoints(
    shapeGrid: Any,
    termPts: Any,
    winRad: Any,
    stepDist: Any,
    maxGridVal: int = ...,
    maxDistC: float = ...,
    distTol: float = ...,
    symFactor: float = ...,
    verbose: bool = ...,
): ...
def CalculateDirectionsAtPoint(pt: Any, shapeGrid: Any, winRad: Any) -> None: ...
def AssignMolFeatsToPoints(
    pts: Any, mol: Any, featFactory: Any, winRad: Any
) -> None: ...
