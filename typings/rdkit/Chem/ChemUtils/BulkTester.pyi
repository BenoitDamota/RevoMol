from rdkit import Chem as Chem
from rdkit.Chem import Randomize as Randomize
from typing import Any

def TestMolecule(mol: Any): ...
def TestSupplier(
    suppl: Any,
    stopAfter: int = ...,
    reportInterval: int = ...,
    reportTo: Any = ...,
    nameProp: str = ...,
): ...
