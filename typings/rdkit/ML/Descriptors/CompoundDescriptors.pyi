from rdkit import RDConfig as RDConfig
from rdkit.ML.Descriptors import Descriptors as Descriptors, Parser as Parser
from rdkit.utils import chemutils as chemutils
from typing import Any, Optional

countOptions: Any

def GetAllDescriptorNames(
    db: Any, tbl1: Any, tbl2: Any, user: str = ..., password: str = ...
): ...

class CompoundDescriptorCalculator(Descriptors.DescriptorCalculator):
    def SUM(self, desc: Any, compos: Any): ...
    def MEAN(self, desc: Any, compos: Any): ...
    def DEV(self, desc: Any, compos: Any): ...
    def MIN(self, desc: Any, compos: Any): ...
    def MAX(self, desc: Any, compos: Any): ...
    nonZeroDescriptors: Any = ...
    requiredDescriptors: Any = ...
    def ProcessSimpleList(self): ...
    def ProcessCompoundList(self) -> None: ...
    atomDict: Any = ...
    def BuildAtomDict(self) -> None: ...
    def CalcSimpleDescriptorsForComposition(
        self, compos: str = ..., composList: Optional[Any] = ...
    ): ...
    def CalcCompoundDescriptorsForComposition(
        self, compos: str = ..., composList: Optional[Any] = ..., propDict: Any = ...
    ): ...
    def CalcDescriptorsForComposition(self, composVect: Any, propDict: Any): ...
    CalcDescriptors: Any = ...
    descriptorNames: Any = ...
    def GetDescriptorNames(self): ...
    simpleList: Any = ...
    compoundList: Any = ...
    dbName: Any = ...
    dbTable: Any = ...
    dbUser: Any = ...
    dbPassword: Any = ...
    def __init__(
        self,
        simpleList: Any,
        compoundList: Optional[Any] = ...,
        dbName: Optional[Any] = ...,
        dbTable: str = ...,
        dbUser: str = ...,
        dbPassword: str = ...,
    ) -> None: ...
