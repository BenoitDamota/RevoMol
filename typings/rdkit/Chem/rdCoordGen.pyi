import Boost.Python

def AddCoords(RDKit) -> Any: ...
def SetDefaultTemplateFileDir(*args, **kwargs) -> Any: ...

class CoordGenParams(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def SetCoordMap(RDKit, boost) -> Any: ...
    @classmethod
    def SetTemplateMol(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
    @property
    def coordgenScaling(self) -> Any: ...
    @coordgenScaling.setter
    def coordgenScaling(self, val: Any) -> None: ...
    @property
    def dbg_useConstrained(self) -> Any: ...
    @dbg_useConstrained.setter
    def dbg_useConstrained(self, val: Any) -> None: ...
    @property
    def dbg_useFixed(self) -> Any: ...
    @dbg_useFixed.setter
    def dbg_useFixed(self, val: Any) -> None: ...
    @property
    def minimizeOnly(self) -> Any: ...
    @minimizeOnly.setter
    def minimizeOnly(self, val: Any) -> None: ...
    @property
    def minimizerPrecision(self) -> Any: ...
    @minimizerPrecision.setter
    def minimizerPrecision(self, val: Any) -> None: ...
    @property
    def sketcherBestPrecision(self) -> Any: ...
    @property
    def sketcherCoarsePrecision(self) -> Any: ...
    @property
    def sketcherQuickPrecision(self) -> Any: ...
    @property
    def sketcherStandardPrecision(self) -> Any: ...
    @property
    def templateFileDir(self) -> Any: ...
    @templateFileDir.setter
    def templateFileDir(self, val: Any) -> None: ...
