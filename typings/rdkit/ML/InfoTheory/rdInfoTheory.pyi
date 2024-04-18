import Boost.Python

BIASCHISQUARE: Any
BIASENTROPY: Any
CHISQUARE: Any
ENTROPY: Any

def ChiSquare(boost) -> Any: ...
def InfoEntropy(boost) -> Any: ...
def InfoGain(boost) -> Any: ...

class BitCorrMatGenerator(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def CollectVotes(RDInfoTheory, boost) -> Any: ...
    @classmethod
    def GetCorrMatrix(RDInfoTheory) -> Any: ...
    @classmethod
    def SetBitList(RDInfoTheory, boost) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class InfoBitRanker(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AccumulateVotes(RDInfoTheory, boost, int) -> Any: ...
    @classmethod
    def GetTopN(RDInfoTheory, int) -> Any: ...
    @classmethod
    def SetBiasList(RDInfoTheory, boost) -> Any: ...
    @classmethod
    def SetMaskBits(RDInfoTheory, boost) -> Any: ...
    @classmethod
    def Tester(RDInfoTheory, boost) -> Any: ...
    @classmethod
    def WriteTopBitsToFile(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class InfoType(Boost.Python.enum):
    BIASCHISQUARE: Any = ...
    BIASENTROPY: Any = ...
    CHISQUARE: Any = ...
    ENTROPY: Any = ...
    names: Any = ...
    values: Any = ...
    __slots__: Any = ...
