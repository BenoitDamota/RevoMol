from rdkit.DataStructs import BitUtils as BitUtils
from rdkit.ML.Data import MLData as MLData
from rdkit.utils import fileutils as fileutils
from typing import Any, Optional

def permutation(nToDo: Any): ...
def WriteData(outFile: Any, varNames: Any, qBounds: Any, examples: Any) -> None: ...
def ReadVars(inFile: Any): ...
def ReadQuantExamples(inFile: Any): ...
def ReadGeneralExamples(inFile: Any): ...
def BuildQuantDataSet(fileName: Any): ...
def BuildDataSet(fileName: Any): ...
def CalcNPossibleUsingMap(
    data: Any,
    order: Any,
    qBounds: Any,
    nQBounds: Optional[Any] = ...,
    silent: bool = ...,
): ...
def WritePickledData(outName: Any, data: Any) -> None: ...
def TakeEnsemble(vect: Any, ensembleIds: Any, isDataVect: bool = ...): ...
def DBToData(
    dbName: Any,
    tableName: Any,
    user: str = ...,
    password: str = ...,
    dupCol: int = ...,
    what: str = ...,
    where: str = ...,
    join: str = ...,
    pickleCol: int = ...,
    pickleClass: Optional[Any] = ...,
    ensembleIds: Optional[Any] = ...,
): ...
def TextToData(reader: Any, ignoreCols: Any = ..., onlyCols: Optional[Any] = ...): ...
def TextFileToData(fName: Any, onlyCols: Optional[Any] = ...): ...
def InitRandomNumbers(seed: Any) -> None: ...
def FilterData(
    inData: Any,
    val: Any,
    frac: Any,
    col: int = ...,
    indicesToUse: Optional[Any] = ...,
    indicesOnly: int = ...,
): ...
def CountResults(inData: Any, col: int = ..., bounds: Optional[Any] = ...): ...
def RandomizeActivities(
    dataSet: Any, shuffle: int = ..., runDetails: Optional[Any] = ...
) -> None: ...