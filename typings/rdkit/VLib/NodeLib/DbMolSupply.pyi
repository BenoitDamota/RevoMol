from rdkit import Chem as Chem, RDConfig as RDConfig
from rdkit.Chem.Suppliers import DbMolSupplier as DbMolSupplier
from rdkit.VLib.Supply import SupplyNode as SupplyNode
from typing import Any

class DbMolSupplyNode(SupplyNode):
    def __init__(self, dbResults: Any, **kwargs: Any) -> None: ...
    def reset(self) -> None: ...
    def next(self): ...

def GetNode(dbName: Any, tableName: Any): ...
