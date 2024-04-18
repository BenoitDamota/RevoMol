from rdkit import Chem as Chem, RDConfig as RDConfig
from rdkit.Chem import AllChem as AllChem, rdChemReactions as rdChemReactions
from typing import Any, Optional

def PreprocessReaction(
    reaction: Any, funcGroupFilename: Optional[Any] = ..., propName: str = ...
): ...
def EnumerateReaction(
    reaction: Any,
    bbLists: Any,
    uniqueProductsOnly: bool = ...,
    funcGroupFilename: Any = ...,
    propName: str = ...,
): ...
