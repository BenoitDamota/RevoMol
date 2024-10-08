from typing import overload
import Boost.Python

def AddMetadataToPNGFile(*args, **kwargs) -> Any: ...
def AddMetadataToPNGString(*args, **kwargs) -> Any: ...
def AtomFromSmarts(*args, **kwargs) -> Any: ...
def AtomFromSmiles(*args, **kwargs) -> Any: ...
def BondFromSmarts(*args, **kwargs) -> Any: ...
def BondFromSmiles(*args, **kwargs) -> Any: ...
@overload
def CanonicalRankAtoms(mol, breakTies = ...) -> Any: ...
@overload
def CanonicalRankAtoms(RDKit) -> Any: ...
@overload
def CanonicalRankAtomsInFragment(mol, atomsToUse = ..., breakTies = ...) -> Any: ...
@overload
def CanonicalRankAtomsInFragment(mol, atomsToUse = ..., breakTies = ...) -> Any: ...
@overload
def CanonicalRankAtomsInFragment(RDKit, boost) -> Any: ...
def CreateAtomBoolPropertyList(*args, **kwargs) -> Any: ...
def CreateAtomDoublePropertyList(*args, **kwargs) -> Any: ...
def CreateAtomIntPropertyList(*args, **kwargs) -> Any: ...
def CreateAtomStringPropertyList(*args, **kwargs) -> Any: ...
def MetadataFromPNGFile(boost) -> Any: ...
def MetadataFromPNGString(boost) -> Any: ...
def MolFragmentToCXSmiles(RDKit, boost) -> Any: ...
def MolFragmentToSmarts(RDKit, boost) -> Any: ...
def MolFragmentToSmiles(RDKit, boost) -> Any: ...
def MolFromFASTA(*args, **kwargs) -> Any: ...
def MolFromHELM(boost) -> Any: ...
def MolFromMol2Block(*args, **kwargs) -> Any: ...
def MolFromMol2File(*args, **kwargs) -> Any: ...
@overload
def MolFromMolBlock(boost) -> Any: ...
@overload
def MolFromMolBlock(boost) -> Any: ...
def MolFromMolFile(*args, **kwargs) -> Any: ...
def MolFromPDBBlock(boost) -> Any: ...
def MolFromPDBFile(*args, **kwargs) -> Any: ...
def MolFromPNGFile(*args, **kwargs) -> Any: ...
def MolFromPNGString(boost) -> Any: ...
def MolFromRDKitSVG(boost) -> Any: ...
def MolFromSequence(boost) -> Any: ...
def MolFromSmarts(boost) -> Any: ...
def MolFromSmiles(*args, **kwargs) -> Any: ...
def MolFromTPLBlock(boost) -> Any: ...
def MolFromTPLFile(*args, **kwargs) -> Any: ...
def MolMetadataToPNGFile(RDKit, boost) -> Any: ...
def MolMetadataToPNGString(RDKit, boost) -> Any: ...
def MolToCXSmiles(RDKit) -> Any: ...
def MolToFASTA(RDKit) -> Any: ...
def MolToHELM(RDKit) -> Any: ...
def MolToMolBlock(RDKit) -> Any: ...
def MolToMolFile(*args, **kwargs) -> Any: ...
def MolToPDBBlock(RDKit) -> Any: ...
def MolToPDBFile(*args, **kwargs) -> Any: ...
def MolToRandomSmilesVect(*args, **kwargs) -> Any: ...
def MolToSequence(RDKit) -> Any: ...
def MolToSmarts(RDKit) -> Any: ...
def MolToSmiles(RDKit) -> Any: ...
def MolToTPLBlock(RDKit) -> Any: ...
def MolToTPLFile(*args, **kwargs) -> Any: ...
def MolToV3KMolBlock(RDKit) -> Any: ...
def MolToV3KMolFile(*args, **kwargs) -> Any: ...
def MolToXYZBlock(RDKit) -> Any: ...
def MolToXYZFile(*args, **kwargs) -> Any: ...
def MolsFromPNGFile(*args, **kwargs) -> Any: ...
def MolsFromPNGString(boost) -> Any: ...
def SmilesMolSupplierFromText(*args, **kwargs) -> Any: ...

class ForwardSDMolSupplier(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetEOFHitOnRead((anonymousnamespace)) -> Any: ...
    @classmethod
    def GetProcessPropertyLists((anonymousnamespace)) -> Any: ...
    @classmethod
    def SetProcessPropertyLists((anonymousnamespace), bool) -> Any: ...
    @classmethod
    def atEnd((anonymousnamespace)) -> Any: ...
    @classmethod
    def __iter__((anonymousnamespace)) -> Any: ...
    @classmethod
    def __next__((anonymousnamespace)) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class MaeMolSupplier(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def atEnd((anonymousnamespace)) -> Any: ...
    @classmethod
    def __iter__((anonymousnamespace)) -> Any: ...
    @classmethod
    def __next__((anonymousnamespace)) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class MultithreadedSDMolSupplier(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetLastItemText(RDKit) -> Any: ...
    @classmethod
    def GetLastRecordId(RDKit) -> Any: ...
    @classmethod
    def GetProcessPropertyLists(RDKit) -> Any: ...
    @classmethod
    def SetProcessPropertyLists(RDKit, bool) -> Any: ...
    @classmethod
    def atEnd(RDKit) -> Any: ...
    @classmethod
    def __iter__(RDKit) -> Any: ...
    @classmethod
    def __next__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class MultithreadedSmilesMolSupplier(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetLastItemText(RDKit) -> Any: ...
    @classmethod
    def GetLastRecordId(RDKit) -> Any: ...
    @classmethod
    def atEnd(RDKit) -> Any: ...
    @classmethod
    def __iter__(RDKit) -> Any: ...
    @classmethod
    def __next__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class PDBWriter(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def NumMols(RDKit) -> Any: ...
    @classmethod
    def close(RDKit) -> Any: ...
    @classmethod
    def flush(RDKit) -> Any: ...
    @classmethod
    def write(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class SDMolSupplier(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetItemText(RDKit, unsignedint) -> Any: ...
    @classmethod
    def GetProcessPropertyLists(RDKit) -> Any: ...
    @classmethod
    def SetData(self, *args, **kwargs) -> Any: ...
    @classmethod
    def SetProcessPropertyLists(RDKit, bool) -> Any: ...
    @classmethod
    def _SetStreamIndices(RDKit, boost) -> Any: ...
    @classmethod
    def atEnd(RDKit) -> Any: ...
    @classmethod
    def reset(RDKit) -> Any: ...
    @classmethod
    def __getitem__(RDKit, int) -> Any: ...
    @classmethod
    def __iter__(RDKit) -> Any: ...
    @classmethod
    def __len__(RDKit) -> Any: ...
    @classmethod
    def __next__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class SDWriter(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetForceV3000(RDKit) -> Any: ...
    @classmethod
    def GetKekulize(RDKit) -> Any: ...
    def GetText(self, *args, **kwargs) -> Any: ...
    @classmethod
    def NumMols(RDKit) -> Any: ...
    @classmethod
    def SetForceV3000(RDKit, bool) -> Any: ...
    @classmethod
    def SetKekulize(RDKit, bool) -> Any: ...
    @classmethod
    def SetProps(RDKit, boost) -> Any: ...
    @classmethod
    def close(RDKit) -> Any: ...
    @classmethod
    def flush(RDKit) -> Any: ...
    @classmethod
    def write(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class SmilesMolSupplier(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetItemText(RDKit, unsignedint) -> Any: ...
    @classmethod
    def SetData(self, *args, **kwargs) -> Any: ...
    @classmethod
    def reset(RDKit) -> Any: ...
    @classmethod
    def __getitem__(RDKit, int) -> Any: ...
    @classmethod
    def __iter__(RDKit) -> Any: ...
    @classmethod
    def __len__(RDKit) -> Any: ...
    @classmethod
    def __next__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class SmilesParserParams(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(self) -> Any: ...
    @property
    def allowCXSMILES(self) -> Any: ...
    @allowCXSMILES.setter
    def allowCXSMILES(self, val: Any) -> None: ...
    @property
    def maxIterations(self) -> Any: ...
    @maxIterations.setter
    def maxIterations(self, val: Any) -> None: ...
    @property
    def parseName(self) -> Any: ...
    @parseName.setter
    def parseName(self, val: Any) -> None: ...
    @property
    def removeHs(self) -> Any: ...
    @removeHs.setter
    def removeHs(self, val: Any) -> None: ...
    @property
    def sanitize(self) -> Any: ...
    @sanitize.setter
    def sanitize(self, val: Any) -> None: ...
    @property
    def strictCXSMILES(self) -> Any: ...
    @strictCXSMILES.setter
    def strictCXSMILES(self, val: Any) -> None: ...
    @property
    def useLegacyStereo(self) -> Any: ...
    @useLegacyStereo.setter
    def useLegacyStereo(self, val: Any) -> None: ...

class SmilesWriter(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def NumMols(RDKit) -> Any: ...
    @classmethod
    def SetProps(RDKit, boost) -> Any: ...
    @classmethod
    def close(RDKit) -> Any: ...
    @classmethod
    def flush(RDKit) -> Any: ...
    @classmethod
    def write(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class TDTMolSupplier(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetItemText(RDKit, unsignedint) -> Any: ...
    @classmethod
    def SetData(self, *args, **kwargs) -> Any: ...
    @classmethod
    def reset(RDKit) -> Any: ...
    @classmethod
    def __getitem__(RDKit, int) -> Any: ...
    @classmethod
    def __iter__(RDKit) -> Any: ...
    @classmethod
    def __len__(RDKit) -> Any: ...
    @classmethod
    def __next__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class TDTWriter(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetNumDigits(RDKit) -> Any: ...
    @classmethod
    def GetWrite2D(RDKit) -> Any: ...
    @classmethod
    def GetWriteNames(RDKit) -> Any: ...
    @classmethod
    def NumMols(RDKit) -> Any: ...
    @classmethod
    def SetNumDigits(RDKit, unsignedint) -> Any: ...
    @classmethod
    def SetProps(RDKit, boost) -> Any: ...
    @classmethod
    def SetWrite2D(RDKit) -> Any: ...
    @classmethod
    def SetWriteNames(RDKit) -> Any: ...
    @classmethod
    def close(RDKit) -> Any: ...
    @classmethod
    def flush(RDKit) -> Any: ...
    @classmethod
    def write(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
