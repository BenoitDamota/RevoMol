from rdkit import Chem as Chem, DataStructs as DataStructs, Geometry as Geometry
from rdkit.Chem import Draw as Draw, rdDepictor as rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D as rdMolDraw2D
from typing import Any, Optional

def GetAtomicWeightsForFingerprint(
    refMol: Any, probeMol: Any, fpFunction: Any, metric: Any = ...
): ...
def GetAtomicWeightsForModel(
    probeMol: Any, fpFunction: Any, predictionFunction: Any
): ...
def GetStandardizedWeights(weights: Any): ...
def GetSimilarityMapFromWeights(
    mol: Any,
    weights: Any,
    colorMap: Optional[Any] = ...,
    scale: int = ...,
    size: Any = ...,
    sigma: Optional[Any] = ...,
    coordScale: float = ...,
    step: float = ...,
    colors: str = ...,
    contourLines: int = ...,
    alpha: float = ...,
    draw2d: Optional[Any] = ...,
    **kwargs: Any,
): ...
def GetSimilarityMapForFingerprint(
    refMol: Any, probeMol: Any, fpFunction: Any, metric: Any = ..., **kwargs: Any
): ...
def GetSimilarityMapForModel(
    probeMol: Any, fpFunction: Any, predictionFunction: Any, **kwargs: Any
): ...

apDict: Any

def GetAPFingerprint(
    mol: Any,
    atomId: int = ...,
    fpType: str = ...,
    nBits: int = ...,
    minLength: int = ...,
    maxLength: int = ...,
    nBitsPerEntry: int = ...,
    **kwargs: Any,
): ...

ttDict: Any

def GetTTFingerprint(
    mol: Any,
    atomId: int = ...,
    fpType: str = ...,
    nBits: int = ...,
    targetSize: int = ...,
    nBitsPerEntry: int = ...,
    **kwargs: Any,
): ...
def GetMorganFingerprint(
    mol: Any,
    atomId: int = ...,
    radius: int = ...,
    fpType: str = ...,
    nBits: int = ...,
    useFeatures: bool = ...,
    **kwargs: Any,
): ...
def GetRDKFingerprint(
    mol: Any,
    atomId: int = ...,
    fpType: str = ...,
    nBits: int = ...,
    minPath: int = ...,
    maxPath: int = ...,
    nBitsPerHash: int = ...,
    **kwargs: Any,
): ...
