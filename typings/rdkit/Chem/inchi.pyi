from typing import Any, Optional

INCHI_AVAILABLE: bool

class InchiReadWriteError(Exception): ...

def MolFromInchi(
    inchi: Any,
    sanitize: bool = ...,
    removeHs: bool = ...,
    logLevel: Optional[Any] = ...,
    treatWarningAsError: bool = ...,
): ...
def MolToInchiAndAuxInfo(
    mol: Any,
    options: str = ...,
    logLevel: Optional[Any] = ...,
    treatWarningAsError: bool = ...,
): ...
def MolBlockToInchiAndAuxInfo(
    molblock: Any,
    options: str = ...,
    logLevel: Optional[Any] = ...,
    treatWarningAsError: bool = ...,
): ...
def MolToInchi(
    mol: Any,
    options: str = ...,
    logLevel: Optional[Any] = ...,
    treatWarningAsError: bool = ...,
): ...
def MolBlockToInchi(
    molblock: Any,
    options: str = ...,
    logLevel: Optional[Any] = ...,
    treatWarningAsError: bool = ...,
): ...
def InchiToInchiKey(inchi: Any): ...
def MolToInchiKey(mol: Any, options: str = ...): ...
