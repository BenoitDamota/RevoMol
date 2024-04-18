from rdkit import Chem as Chem
from rdkit.Chem import (
    ChemicalFeatures as ChemicalFeatures,
    ChemicalForceFields as ChemicalForceFields,
)
from rdkit.Chem.Pharm3D import ExcludedVolume as ExcludedVolume
from rdkit.ML.Data import Stats as Stats
from typing import Any, Optional

logger: Any
defaultFeatLength: float

def GetAtomHeavyNeighbors(atom: Any): ...
def ReplaceGroup(
    match: Any,
    bounds: Any,
    slop: float = ...,
    useDirs: bool = ...,
    dirLength: Any = ...,
): ...
def EmbedMol(
    mol: Any,
    bm: Any,
    atomMatch: Optional[Any] = ...,
    weight: float = ...,
    randomSeed: int = ...,
    excludedVolumes: Optional[Any] = ...,
) -> None: ...
def AddExcludedVolumes(bm: Any, excludedVolumes: Any, smoothIt: bool = ...): ...
def UpdatePharmacophoreBounds(
    bm: Any,
    atomMatch: Any,
    pcophore: Any,
    useDirs: bool = ...,
    dirLength: Any = ...,
    mol: Optional[Any] = ...,
): ...
def EmbedPharmacophore(
    mol: Any,
    atomMatch: Any,
    pcophore: Any,
    randomSeed: int = ...,
    count: int = ...,
    smoothFirst: bool = ...,
    silent: bool = ...,
    bounds: Optional[Any] = ...,
    excludedVolumes: Optional[Any] = ...,
    targetNumber: int = ...,
    useDirs: bool = ...,
): ...
def isNaN(v: Any): ...
def OptimizeMol(
    mol: Any,
    bm: Any,
    atomMatches: Optional[Any] = ...,
    excludedVolumes: Optional[Any] = ...,
    forceConstant: float = ...,
    maxPasses: int = ...,
    verbose: bool = ...,
): ...
def EmbedOne(
    mol: Any,
    name: Any,
    match: Any,
    pcophore: Any,
    count: int = ...,
    silent: int = ...,
    **kwargs: Any,
): ...
def MatchPharmacophoreToMol(mol: Any, featFactory: Any, pcophore: Any): ...
def MatchFeatsToMol(mol: Any, featFactory: Any, features: Any): ...
def CombiEnum(sequence: Any) -> None: ...
def DownsampleBoundsMatrix(bm: Any, indices: Any, maxThresh: float = ...): ...
def CoarseScreenPharmacophore(
    atomMatch: Any, bounds: Any, pcophore: Any, verbose: bool = ...
): ...
def Check2DBounds(atomMatch: Any, mol: Any, pcophore: Any): ...
def ConstrainedEnum(
    matches: Any,
    mol: Any,
    pcophore: Any,
    bounds: Any,
    use2DLimits: bool = ...,
    index: int = ...,
    soFar: Any = ...,
) -> None: ...
def MatchPharmacophore(
    matches: Any,
    bounds: Any,
    pcophore: Any,
    useDownsampling: bool = ...,
    use2DLimits: bool = ...,
    mol: Optional[Any] = ...,
    excludedVolumes: Optional[Any] = ...,
    useDirs: bool = ...,
): ...
def GetAllPharmacophoreMatches(
    matches: Any,
    bounds: Any,
    pcophore: Any,
    useDownsampling: int = ...,
    progressCallback: Optional[Any] = ...,
    use2DLimits: bool = ...,
    mol: Optional[Any] = ...,
    verbose: bool = ...,
): ...
def ComputeChiralVolume(mol: Any, centerIdx: Any, confId: int = ...): ...
