from rdkit import Chem as Chem
from typing import Any

log: Any

class MetalDisconnector:
    def __init__(self) -> None: ...
    def __call__(self, mol: Any): ...
    def disconnect(self, mol: Any): ...
