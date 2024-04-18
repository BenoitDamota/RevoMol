from typing import overload
import Boost.Python

SANITIZE_ADJUST_REACTANTS: Any
SANITIZE_ALL: Any
SANITIZE_ATOM_MAPS: Any
SANITIZE_MERGEHS: Any
SANITIZE_NONE: Any
SANITIZE_RGROUP_NAMES: Any

def Compute2DCoordsForReaction(RDKit) -> Any: ...
def CreateDifferenceFingerprintForReaction(RDKit) -> Any: ...
def CreateStructuralFingerprintForReaction(RDKit) -> Any: ...
def EnumerateLibraryCanSerialize(*args, **kwargs) -> Any: ...
def GetChemDrawRxnAdjustParams(*args, **kwargs) -> Any: ...
def GetDefaultAdjustParams(*args, **kwargs) -> Any: ...
def HasAgentTemplateSubstructMatch(*args, **kwargs) -> Any: ...
def HasProductTemplateSubstructMatch(*args, **kwargs) -> Any: ...
def HasReactantTemplateSubstructMatch(*args, **kwargs) -> Any: ...
def HasReactionAtomMapping(RDKit) -> Any: ...
def HasReactionSubstructMatch(*args, **kwargs) -> Any: ...
def IsReactionTemplateMoleculeAgent(RDKit, double) -> Any: ...
def MatchOnlyAtRgroupsAdjustParams(*args, **kwargs) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(RDKit) -> Any: ...
def ReactionFromMolecule(RDKit) -> Any: ...
def ReactionFromPNGFile(*args, **kwargs) -> Any: ...
def ReactionFromPNGString(*args, **kwargs) -> Any: ...
def ReactionFromRxnBlock(*args, **kwargs) -> Any: ...
def ReactionFromRxnFile(*args, **kwargs) -> Any: ...
def ReactionFromSmarts(*args, **kwargs) -> Any: ...
def ReactionMetadataToPNGFile(RDKit, boost) -> Any: ...
def ReactionMetadataToPNGString(RDKit, boost) -> Any: ...
def ReactionToMolecule(RDKit) -> Any: ...
def ReactionToRxnBlock(RDKit) -> Any: ...
def ReactionToSmarts(RDKit) -> Any: ...
def ReactionToSmiles(RDKit) -> Any: ...
def ReduceProductToSideChains(boost) -> Any: ...
def RemoveMappingNumbersFromReactions(RDKit) -> Any: ...
def SanitizeRxn(*args, **kwargs) -> Any: ...
def UpdateProductsStereochemistry(RDKit) -> Any: ...

class CartesianProductStrategy(EnumerationStrategyBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class ChemicalReaction(Boost.Python.instance):
    __instance_size__: Any = ...
    __safe_for_unpickling__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AddAgentTemplate(RDKit, boost) -> Any: ...
    @classmethod
    def AddProductTemplate(RDKit, boost) -> Any: ...
    @classmethod
    def AddReactantTemplate(RDKit, boost) -> Any: ...
    @classmethod
    def AddRecursiveQueriesToReaction(RDKit) -> Any: ...
    @classmethod
    def ClearComputedProps(RDKit) -> Any: ...
    @classmethod
    def ClearProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetAgentTemplate(RDKit, unsignedint) -> Any: ...
    @classmethod
    def GetAgents(RDKit) -> Any: ...
    @classmethod
    def GetBoolProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetDoubleProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetIntProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetNumAgentTemplates(RDKit) -> Any: ...
    @classmethod
    def GetNumProductTemplates(RDKit) -> Any: ...
    @classmethod
    def GetNumReactantTemplates(RDKit) -> Any: ...
    @classmethod
    def GetProductTemplate(RDKit, unsignedint) -> Any: ...
    @classmethod
    def GetProducts(RDKit) -> Any: ...
    @classmethod
    def GetProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetPropNames(RDKit) -> Any: ...
    @classmethod
    def GetPropsAsDict(RDKit) -> Any: ...
    @classmethod
    def GetReactantTemplate(RDKit, unsignedint) -> Any: ...
    @classmethod
    def GetReactants(RDKit) -> Any: ...
    @classmethod
    def GetReactingAtoms(RDKit) -> Any: ...
    @classmethod
    def GetUnsignedProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def HasProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def Initialize(RDKit) -> Any: ...
    @classmethod
    def IsInitialized(RDKit) -> Any: ...
    @classmethod
    def IsMoleculeAgent(self, *args, **kwargs) -> Any: ...
    @classmethod
    def IsMoleculeProduct(self, *args, **kwargs) -> Any: ...
    @classmethod
    def IsMoleculeReactant(self, *args, **kwargs) -> Any: ...
    @classmethod
    def RemoveAgentTemplates(RDKit) -> Any: ...
    @classmethod
    def RemoveUnmappedProductTemplates(RDKit) -> Any: ...
    @classmethod
    def RemoveUnmappedReactantTemplates(RDKit) -> Any: ...
    @classmethod
    def RunReactant(RDKit, boost, unsignedint) -> Any: ...
    @classmethod
    @overload
    def RunReactants(RDKit, boost) -> Any: ...
    @overload
    def RunReactants(RDKit, boost) -> Any: ...
    @classmethod
    def SetBoolProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def SetDoubleProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def SetIntProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def SetProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    def SetUnsignedProp(self, *args, **kwargs) -> Any: ...
    @classmethod
    @overload
    def ToBinary(RDKit) -> Any: ...
    @overload
    def ToBinary(RDKit, unsignedint) -> Any: ...
    @classmethod
    def Validate(RDKit) -> Any: ...
    @classmethod
    def _getImplicitPropertiesFlag(RDKit) -> Any: ...
    @classmethod
    def _setImplicitPropertiesFlag(RDKit, bool) -> Any: ...
    @classmethod
    def __getinitargs__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class EnumerateLibrary(EnumerateLibraryBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetReagents(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class EnumerateLibraryBase(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetEnumerator(RDKit) -> Any: ...
    @classmethod
    def GetPosition(RDKit) -> Any: ...
    @classmethod
    def GetReaction(RDKit) -> Any: ...
    @classmethod
    def GetState(RDKit) -> Any: ...
    @classmethod
    def InitFromString(self, *args, **kwargs) -> Any: ...
    @classmethod
    def ResetState(RDKit) -> Any: ...
    @classmethod
    def Serialize(RDKit) -> Any: ...
    @classmethod
    def SetState(self, *args, **kwargs) -> Any: ...
    @classmethod
    def next(RDKit) -> Any: ...
    @classmethod
    def nextSmiles(RDKit) -> Any: ...
    @classmethod
    def __bool__(RDKit) -> Any: ...
    @classmethod
    def __iter__(boost) -> Any: ...
    @classmethod
    def __next__(RDKit) -> Any: ...
    @classmethod
    def __nonzero__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class EnumerationParams(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(self) -> Any: ...
    @property
    def reagentMaxMatchCount(self) -> Any: ...
    @reagentMaxMatchCount.setter
    def reagentMaxMatchCount(self, val: Any) -> None: ...
    @property
    def sanePartialProducts(self) -> Any: ...
    @sanePartialProducts.setter
    def sanePartialProducts(self, val: Any) -> None: ...

class EnumerationStrategyBase(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetNumPermutations(RDKit) -> Any: ...
    @classmethod
    def GetPosition(RDKit) -> Any: ...
    @classmethod
    def Initialize(self, *args, **kwargs) -> Any: ...
    @classmethod
    def Skip(RDKit, unsignedlong) -> Any: ...
    @classmethod
    def Type(RDKit) -> Any: ...
    @classmethod
    @overload
    def next(RDKit) -> Any: ...
    @overload
    def next(RDKit) -> Any: ...
    @classmethod
    def __bool__(RDKit) -> Any: ...
    @classmethod
    @overload
    def __copy__(RDKit) -> Any: ...
    @overload
    def __copy__(RDKit) -> Any: ...
    @classmethod
    @overload
    def __next__(RDKit) -> Any: ...
    @overload
    def __next__(RDKit) -> Any: ...
    @classmethod
    def __nonzero__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class EvenSamplePairsStrategy(EnumerationStrategyBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def Stats(RDKit) -> Any: ...
    @classmethod
    def __copy__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class FingerprintType(Boost.Python.enum):
    AtomPairFP: Any = ...
    MorganFP: Any = ...
    PatternFP: Any = ...
    RDKitFP: Any = ...
    TopologicalTorsion: Any = ...
    names: Any = ...
    values: Any = ...
    __slots__: Any = ...

class MOL_SPTR_VECT(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def append(self, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __contains__(self, other) -> Any: ...
    @classmethod
    def __delitem__(self, other) -> Any: ...
    @classmethod
    def __getitem__(self, index) -> Any: ...
    @classmethod
    def __iter__(boost, std) -> Any: ...
    @classmethod
    def __len__(self) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
    @classmethod
    def __setitem__(self, index, object) -> Any: ...

class RandomSampleAllBBsStrategy(EnumerationStrategyBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class RandomSampleStrategy(EnumerationStrategyBase):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(RDKit) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class ReactionFingerprintParams(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(self) -> Any: ...
    @property
    def agentWeight(self) -> Any: ...
    @agentWeight.setter
    def agentWeight(self, val: Any) -> None: ...
    @property
    def bitRatioAgents(self) -> Any: ...
    @bitRatioAgents.setter
    def bitRatioAgents(self, val: Any) -> None: ...
    @property
    def fpSize(self) -> Any: ...
    @fpSize.setter
    def fpSize(self, val: Any) -> None: ...
    @property
    def fpType(self) -> Any: ...
    @fpType.setter
    def fpType(self, val: Any) -> None: ...
    @property
    def includeAgents(self) -> Any: ...
    @includeAgents.setter
    def includeAgents(self, val: Any) -> None: ...
    @property
    def nonAgentWeight(self) -> Any: ...
    @nonAgentWeight.setter
    def nonAgentWeight(self, val: Any) -> None: ...

class SanitizeFlags(Boost.Python.enum):
    SANITIZE_ADJUST_REACTANTS: Any = ...
    SANITIZE_ALL: Any = ...
    SANITIZE_ATOM_MAPS: Any = ...
    SANITIZE_MERGEHS: Any = ...
    SANITIZE_NONE: Any = ...
    SANITIZE_RGROUP_NAMES: Any = ...
    names: Any = ...
    values: Any = ...
    __slots__: Any = ...

class VectMolVect(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def append(self, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __contains__(self, other) -> Any: ...
    @classmethod
    def __delitem__(self, other) -> Any: ...
    @classmethod
    def __getitem__(self, index) -> Any: ...
    @classmethod
    def __iter__(self) -> Any: ...
    @classmethod
    def __len__(self) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
    @classmethod
    def __setitem__(self, index, object) -> Any: ...

class VectSizeT(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def append(self, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __contains__(self, other) -> Any: ...
    @classmethod
    def __delitem__(self, other) -> Any: ...
    @classmethod
    def __getitem__(self, index) -> Any: ...
    @classmethod
    def __iter__(boost, std) -> Any: ...
    @classmethod
    def __len__(self) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
    @classmethod
    def __setitem__(self, index, object) -> Any: ...

class VectorOfStringVectors(Boost.Python.instance):
    __instance_size__: Any = ...
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def append(self, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __contains__(self, other) -> Any: ...
    @classmethod
    def __delitem__(self, other) -> Any: ...
    @classmethod
    def __getitem__(self, index) -> Any: ...
    @classmethod
    def __iter__(self) -> Any: ...
    @classmethod
    def __len__(self) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
    @classmethod
    def __setitem__(self, index, object) -> Any: ...
