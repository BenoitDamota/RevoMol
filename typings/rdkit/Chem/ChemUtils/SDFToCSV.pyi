from rdkit import Chem as Chem
from typing import Any, Optional

def Convert(
    suppl: Any,
    outFile: Any,
    keyCol: Optional[Any] = ...,
    stopAfter: int = ...,
    includeChirality: bool = ...,
    smilesFrom: str = ...,
) -> None: ...
def initParser(): ...
def existingFile(filename: Any): ...
def main() -> None: ...
