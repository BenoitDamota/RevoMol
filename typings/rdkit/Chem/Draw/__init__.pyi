from rdkit.Chem.Draw.rdMolDraw2D import *
from collections import namedtuple
from rdkit import Chem as Chem, rdBase as rdBase
from rdkit.Chem import rdDepictor as rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D as rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import (
    DrawingOptions as DrawingOptions,
    MolDrawing as MolDrawing,
)
from typing import Any, Optional

def MolToImage(
    mol: Any,
    size: Any = ...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    fitImage: bool = ...,
    options: Optional[Any] = ...,
    canvas: Optional[Any] = ...,
    **kwargs: Any,
): ...
def MolToFile(
    mol: Any,
    filename: Any,
    size: Any = ...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    imageType: Optional[Any] = ...,
    fitImage: bool = ...,
    options: Optional[Any] = ...,
    **kwargs: Any,
) -> None: ...
def MolToImageFile(
    mol: Any,
    filename: Any,
    size: Any = ...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    **kwargs: Any,
) -> None: ...
def ShowMol(
    mol: Any,
    size: Any = ...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    title: str = ...,
    stayInFront: bool = ...,
    **kwargs: Any,
) -> None: ...
def MolToMPL(
    mol: Any,
    size: Any = ...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    imageType: Optional[Any] = ...,
    fitImage: bool = ...,
    options: Optional[Any] = ...,
    **kwargs: Any,
): ...
def calcAtomGaussians(
    mol: Any, a: float = ..., step: float = ..., weights: Optional[Any] = ...
): ...
def MolsToImage(
    mols: Any, subImgSize: Any = ..., legends: Optional[Any] = ..., **kwargs: Any
): ...
def MolsToGridImage(
    mols: Any,
    molsPerRow: int = ...,
    subImgSize: Any = ...,
    legends: Optional[Any] = ...,
    highlightAtomLists: Optional[Any] = ...,
    highlightBondLists: Optional[Any] = ...,
    useSVG: bool = ...,
    returnPNG: bool = ...,
    **kwargs: Any,
): ...
def ReactionToImage(
    rxn: Any,
    subImgSize: Any = ...,
    useSVG: bool = ...,
    drawOptions: Optional[Any] = ...,
    returnPNG: bool = ...,
    **kwargs: Any,
): ...
def MolToQPixmap(
    mol: Any,
    size: Any = ...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    fitImage: bool = ...,
    options: Optional[Any] = ...,
    **kwargs: Any,
): ...
def DrawMorganBit(
    mol: Any, bitId: Any, bitInfo: Any, whichExample: int = ..., **kwargs: Any
): ...
def DrawMorganBits(tpls: Any, **kwargs: Any): ...

FingerprintEnv = namedtuple(
    "FingerprintEnv",
    [
        "submol",
        "highlightAtoms",
        "atomColors",
        "highlightBonds",
        "bondColors",
        "highlightRadii",
    ],
)

def DrawMorganEnvs(
    envs: Any,
    molsPerRow: int = ...,
    subImgSize: Any = ...,
    baseRad: float = ...,
    useSVG: bool = ...,
    aromaticColor: Any = ...,
    ringColor: Any = ...,
    centerColor: Any = ...,
    extraColor: Any = ...,
    legends: Optional[Any] = ...,
    drawOptions: Optional[Any] = ...,
    **kwargs: Any,
): ...
def DrawMorganEnv(
    mol: Any,
    atomId: Any,
    radius: Any,
    molSize: Any = ...,
    baseRad: float = ...,
    useSVG: bool = ...,
    aromaticColor: Any = ...,
    ringColor: Any = ...,
    centerColor: Any = ...,
    extraColor: Any = ...,
    drawOptions: Optional[Any] = ...,
    **kwargs: Any,
): ...
def DrawRDKitBits(tpls: Any, **kwargs: Any): ...
def DrawRDKitBit(
    mol: Any, bitId: Any, bitInfo: Any, whichExample: int = ..., **kwargs: Any
): ...
def DrawRDKitEnvs(
    envs: Any,
    molsPerRow: int = ...,
    subImgSize: Any = ...,
    baseRad: float = ...,
    useSVG: bool = ...,
    aromaticColor: Any = ...,
    extraColor: Any = ...,
    nonAromaticColor: Optional[Any] = ...,
    legends: Optional[Any] = ...,
    drawOptions: Optional[Any] = ...,
    **kwargs: Any,
): ...
def DrawRDKitEnv(
    mol: Any,
    bondPath: Any,
    molSize: Any = ...,
    baseRad: float = ...,
    useSVG: bool = ...,
    aromaticColor: Any = ...,
    extraColor: Any = ...,
    nonAromaticColor: Optional[Any] = ...,
    drawOptions: Optional[Any] = ...,
    **kwargs: Any,
): ...
