from rdkit import RDConfig as RDConfig
from typing import Any, Optional

def SmilesToGif(
    smiles: Any,
    fileNames: Any,
    size: Any = ...,
    cmd: Optional[Any] = ...,
    dblSize: int = ...,
    frame: int = ...,
): ...
def SmilesToImage(smiles: Any, **kwargs: Any): ...
