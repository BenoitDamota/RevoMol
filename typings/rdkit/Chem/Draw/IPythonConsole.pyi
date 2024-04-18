from io import StringIO as StringIO
from rdkit import Chem as Chem
from rdkit.Chem import (
    Draw as Draw,
    rdChemReactions as rdChemReactions,
    rdchem as rdchem,
)
from rdkit.Chem.Draw import rdMolDraw2D as rdMolDraw2D
from typing import Any, Optional

molSize: Any
highlightSubstructs: bool
kekulizeStructures: bool
highlightByReactant: bool
ipython_useSVG: bool
ipython_3d: bool
molSize_3d: Any
drawing_type_3d: str
bgcolor_3d: str
drawOptions: Any

def addMolToView(
    mol: Any, view: Any, confId: int = ..., drawAs: Optional[Any] = ...
) -> None: ...
def drawMol3D(
    m: Any,
    view: Optional[Any] = ...,
    confId: int = ...,
    drawAs: Optional[Any] = ...,
    bgColor: Optional[Any] = ...,
    size: Optional[Any] = ...,
): ...
def display_pil_image(img: Any): ...
def ShowMols(mols: Any, maxMols: int = ..., **kwargs: Any): ...
def DrawMorganBit(
    mol: Any, bitId: Any, bitInfo: Any, drawOptions: Any = ..., **kwargs: Any
): ...
def DrawMorganBits(*args: Any, drawOptions: Any = ..., **kwargs: Any): ...
def DrawRDKitBit(
    mol: Any, bitId: Any, bitInfo: Any, drawOptions: Any = ..., **kwargs: Any
): ...
def DrawRDKitBits(*args: Any, drawOptions: Any = ..., **kwargs: Any): ...
def EnableSubstructMatchRendering() -> None: ...
def InstallIPythonRenderer(): ...
def DisableSubstructMatchRendering() -> None: ...
def UninstallIPythonRenderer() -> None: ...
