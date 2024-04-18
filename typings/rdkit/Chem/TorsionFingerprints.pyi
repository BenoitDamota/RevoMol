from rdkit import (
    Chem as Chem,
    Geometry as Geometry,
    RDConfig as RDConfig,
    rdBase as rdBase,
)
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors, rdchem as rdchem
from typing import Any, Optional

def CalculateTorsionLists(
    mol: Any, maxDev: str = ..., symmRadius: int = ..., ignoreColinearBonds: bool = ...
): ...
def CalculateTorsionAngles(
    mol: Any, tors_list: Any, tors_list_rings: Any, confId: int = ...
): ...
def CalculateTorsionWeights(
    mol: Any, aid1: int = ..., aid2: int = ..., ignoreColinearBonds: bool = ...
): ...
def CalculateTFD(torsions1: Any, torsions2: Any, weights: Optional[Any] = ...): ...
def GetTFDBetweenConformers(
    mol: Any,
    confIds1: Any,
    confIds2: Any,
    useWeights: bool = ...,
    maxDev: str = ...,
    symmRadius: int = ...,
    ignoreColinearBonds: bool = ...,
): ...
def GetTFDBetweenMolecules(
    mol1: Any,
    mol2: Any,
    confId1: int = ...,
    confId2: int = ...,
    useWeights: bool = ...,
    maxDev: str = ...,
    symmRadius: int = ...,
    ignoreColinearBonds: bool = ...,
): ...
def GetTFDMatrix(
    mol: Any,
    useWeights: bool = ...,
    maxDev: str = ...,
    symmRadius: int = ...,
    ignoreColinearBonds: bool = ...,
): ...
