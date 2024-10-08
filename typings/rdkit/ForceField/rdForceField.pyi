import Boost.Python

class ForceField(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def AddDistanceConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def AddExtraPoint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def AddFixedPoint(ForceFields, unsignedint) -> Any: ...
    @classmethod
    def CalcEnergy(ForceFields) -> Any: ...
    @classmethod
    def CalcGrad(ForceFields) -> Any: ...
    @classmethod
    def Dimension(ForceFields) -> Any: ...
    @classmethod
    def GetExtraPointPos(ForceFields, unsignedint) -> Any: ...
    @classmethod
    def Initialize(ForceFields) -> Any: ...
    @classmethod
    def MMFFAddAngleConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def MMFFAddDistanceConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def MMFFAddPositionConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def MMFFAddTorsionConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def Minimize(ForceFields) -> Any: ...
    @classmethod
    def MinimizeTrajectory(self, *args, **kwargs) -> Any: ...
    @classmethod
    def NumPoints(ForceFields) -> Any: ...
    @classmethod
    def Positions(ForceFields) -> Any: ...
    @classmethod
    def UFFAddAngleConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def UFFAddDistanceConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def UFFAddPositionConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def UFFAddTorsionConstraint(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...

class MMFFMolProperties(Boost.Python.instance):
    @classmethod
    def __init__(self, *args, **kwargs) -> None: ...
    @classmethod
    def GetMMFFAngleBendParams(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFAtomType(ForceFields, unsignedint) -> Any: ...
    @classmethod
    def GetMMFFBondStretchParams(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFFormalCharge(ForceFields, unsignedint) -> Any: ...
    @classmethod
    def GetMMFFOopBendParams(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFPartialCharge(ForceFields, unsignedint) -> Any: ...
    @classmethod
    def GetMMFFStretchBendParams(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFTorsionParams(self, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFVdWParams(self, *args, **kwargs) -> Any: ...
    @classmethod
    def SetMMFFAngleTerm(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFBondTerm(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFDielectricConstant(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFDielectricModel(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFEleTerm(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFOopTerm(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFStretchBendTerm(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFTorsionTerm(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFVariant(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFVdWTerm(ForceFields) -> Any: ...
    @classmethod
    def SetMMFFVerbosity(ForceFields) -> Any: ...
    @classmethod
    def __reduce__(self) -> Any: ...
