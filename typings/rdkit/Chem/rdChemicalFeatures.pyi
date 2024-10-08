import Boost.Python

class FreeChemicalFeature(Boost.Python.instance):
    __instance_size__: Any = ...
    __safe_for_unpickling__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetFamily(ChemicalFeatures) -> Any: ...
    @classmethod
    def GetId(ChemicalFeatures) -> Any: ...
    @classmethod
    def GetPos(ChemicalFeatures) -> Any: ...
    @classmethod
    def GetType(ChemicalFeatures) -> Any: ...
    @classmethod
    def SetFamily(self, *args, **kwargs) -> Any: ...
    @classmethod
    def SetId(ChemicalFeatures, int) -> Any: ...
    @classmethod
    def SetPos(ChemicalFeatures, RDGeom) -> Any: ...
    @classmethod
    def SetType(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __getinitargs__(ChemicalFeatures) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
