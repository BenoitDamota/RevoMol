from rdkit.sping.pid import *
from . import psmetrics as psmetrics
from typing import Any, Optional

class PostScriptLevelException(ValueError): ...

linesep: str
PiddleLegalFonts: Any
Roman: str
Bold: str
Italic: str
PSFontMapStdEnc: Any
PSFontMapLatin1Enc: Any

def latin1FontEncoding(fontname: Any): ...
def dashLineDefinition(): ...

class PsDSC:
    def __init__(self) -> None: ...
    def documentHeader(self): ...
    def boundingBoxStr(self, x0: Any, y0: Any, x1: Any, y1: Any): ...
    inPageFlag: int = ...
    def BeginPageStr(self, pageSetupStr: Any, pageName: Optional[Any] = ...): ...
    def EndPageStr(self): ...

class EpsDSC(PsDSC):
    def __init__(self) -> None: ...
    def documentHeader(self): ...

class PSCanvas(Canvas):
    filename: Any = ...
    drawImage: Any = ...
    code: Any = ...
    dsc: Any = ...
    defaultFont: Any = ...
    fontMapEncoding: Any = ...
    pageNum: int = ...
    def __init__(
        self,
        size: Any = ...,
        name: str = ...,
        PostScriptLevel: int = ...,
        fontMapEncoding: Any = ...,
    ) -> None: ...
    def psBeginDocument(self) -> None: ...
    def psEndDocument(self) -> None: ...
    def psBeginPage(self, pageName: Optional[Any] = ...) -> None: ...
    def psEndPage(self) -> None: ...
    def nextPage(self) -> None: ...
    def clear(self) -> None: ...
    def resetToDefaults(self) -> None: ...
    def flush(self) -> None: ...
    def save(self, file: Optional[Any] = ..., format: Optional[Any] = ...) -> None: ...
    def stringWidth(self, s: Any, font: Optional[Any] = ...): ...
    def fontAscent(self, font: Optional[Any] = ...): ...
    def fontDescent(self, font: Optional[Any] = ...): ...
    def drawLine(
        self,
        x1: Any,
        y1: Any,
        x2: Any,
        y2: Any,
        color: Optional[Any] = ...,
        width: Optional[Any] = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawLines(
        self,
        lineList: Any,
        color: Optional[Any] = ...,
        width: Optional[Any] = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawString(
        self,
        s: Any,
        x: Any,
        y: Any,
        font: Optional[Any] = ...,
        color: Optional[Any] = ...,
        angle: int = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawCurve(
        self,
        x1: Any,
        y1: Any,
        x2: Any,
        y2: Any,
        x3: Any,
        y3: Any,
        x4: Any,
        y4: Any,
        edgeColor: Optional[Any] = ...,
        edgeWidth: Optional[Any] = ...,
        fillColor: Optional[Any] = ...,
        closed: int = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawRoundRect(
        self,
        x1: Any,
        y1: Any,
        x2: Any,
        y2: Any,
        rx: int = ...,
        ry: int = ...,
        edgeColor: Optional[Any] = ...,
        edgeWidth: Optional[Any] = ...,
        fillColor: Optional[Any] = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawEllipse(
        self,
        x1: Any,
        y1: Any,
        x2: Any,
        y2: Any,
        edgeColor: Optional[Any] = ...,
        edgeWidth: Optional[Any] = ...,
        fillColor: Optional[Any] = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawArc(
        self,
        x1: Any,
        y1: Any,
        x2: Any,
        y2: Any,
        startAng: int = ...,
        extent: int = ...,
        edgeColor: Optional[Any] = ...,
        edgeWidth: Optional[Any] = ...,
        fillColor: Optional[Any] = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawPolygon(
        self,
        pointlist: Any,
        edgeColor: Optional[Any] = ...,
        edgeWidth: Optional[Any] = ...,
        fillColor: Optional[Any] = ...,
        closed: int = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawFigure(
        self,
        partList: Any,
        edgeColor: Optional[Any] = ...,
        edgeWidth: Optional[Any] = ...,
        fillColor: Optional[Any] = ...,
        closed: int = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
