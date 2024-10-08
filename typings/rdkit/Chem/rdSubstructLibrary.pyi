from typing import overload
import Boost.Python

def AddPatterns(RDKit) -> Any: ...
def SubstructLibraryCanSerialize(*args, **kwargs) -> Any: ...

class CachedMolHolder(MolHolderBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AddBinary(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class CachedSmilesMolHolder(MolHolderBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AddSmiles(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class CachedTrustedSmilesMolHolder(MolHolderBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AddSmiles(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class FPHolderBase(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AddFingerprint(RDKit, ExplicitBitVect) -> Any: ...
    @classmethod
    def AddMol(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetFingerprint(RDKit, unsignedint) -> Any: ...
    @classmethod
    def MakeFingerprint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def PassesFilter(RDKit, unsignedint, ExplicitBitVect) -> Any: ...
    @classmethod
    def __len__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class MolHolder(MolHolderBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class MolHolderBase(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AddMol(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMol(RDKit, unsignedint) -> Any: ...
    @classmethod
    @overload
    def __len__(RDKit) -> Any: ...
    @overload
    def __len__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class PatternHolder(FPHolderBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class SubstructLibrary(Boost.Python.instance):
    __instance_size__: Any = ...
    __safe_for_unpickling__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AddMol(self, *args, **kwargs) -> Any: ...
    @classmethod
    def CountMatches(RDKit) -> Any: ...
    @classmethod
    def GetFpHolder(RDKit) -> Any: ...
    @classmethod
    def GetMatches(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMol(RDKit, unsignedint) -> Any: ...
    @classmethod
    def GetMolHolder(RDKit) -> Any: ...
    @classmethod
    def HasMatch(self, *args, **kwargs) -> Any: ...
    @classmethod
    @overload
    def InitFromStream(stream) -> Any: ...
    @overload
    def InitFromStream(f) -> Any: ...
    @overload
    def InitFromStream(RDKit, boost) -> Any: ...
    @classmethod
    def Serialize(RDKit) -> Any: ...
    @classmethod
    @overload
    def ToStream(stream) -> Any: ...
    @overload
    def ToStream(stream) -> Any: ...
    @overload
    def ToStream(RDKit, boost) -> Any: ...
    @classmethod
    def __getinitargs__(RDKit) -> Any: ...
    @classmethod
    def __len__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
