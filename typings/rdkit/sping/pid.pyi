from rdkit.sping.colors import *
from typing import Any, Optional

__version_maj_number__: str
__version_min_number__: str
inch: int
cm: Any

class StateSaver:
    canvas: Any = ...
    defaultLineColor: Any = ...
    defaultFillColor: Any = ...
    defaultLineWidth: Any = ...
    defaultFont: Any = ...
    def __init__(self, canvas: Any) -> None: ...
    def __del__(self) -> None: ...

class Font:
    def __init__(
        self,
        size: int = ...,
        bold: int = ...,
        italic: int = ...,
        underline: int = ...,
        face: Optional[Any] = ...,
    ) -> None: ...
    def __cmp__(self, other: Any): ...
    def __setattr__(self, name: Any, value: Any) -> None: ...

figureLine: int
figureArc: int
figureCurve: int
keyBksp: str
keyDel: str
keyLeft: str
keyRight: str
keyUp: str
keyDown: str
keyPgUp: str
keyPgDn: str
keyHome: str
keyEnd: str
keyClear: str
keyTab: str
modShift: int
modControl: int

class Canvas:
    defaultLineColor: Any = ...
    defaultFillColor: Any = ...
    defaultLineWidth: int = ...
    defaultFont: Any = ...
    onClick: Any = ...
    onOver: Any = ...
    onKey: Any = ...
    def __init__(self, size: Any = ..., name: str = ...) -> None: ...
    def getSize(self): ...
    def isInteractive(self): ...
    def canUpdate(self): ...
    def clear(self) -> None: ...
    def flush(self) -> None: ...
    def save(self, file: Optional[Any] = ..., format: Optional[Any] = ...) -> None: ...
    def setInfoLine(self, s: Any) -> None: ...
    def stringBox(self, s: Any, font: Optional[Any] = ...): ...
    def stringWidth(self, s: Any, font: Optional[Any] = ...) -> None: ...
    def fontHeight(self, font: Optional[Any] = ...): ...
    def fontAscent(self, font: Optional[Any] = ...) -> None: ...
    def fontDescent(self, font: Optional[Any] = ...) -> None: ...
    def arcPoints(
        self, x1: Any, y1: Any, x2: Any, y2: Any, startAng: int = ..., extent: int = ...
    ): ...
    def curvePoints(
        self, x1: Any, y1: Any, x2: Any, y2: Any, x3: Any, y3: Any, x4: Any, y4: Any
    ): ...
    def drawMultiLineString(
        self,
        s: Any,
        x: Any,
        y: Any,
        font: Optional[Any] = ...,
        color: Optional[Any] = ...,
        angle: int = ...,
        **kwargs: Any,
    ) -> None: ...
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
    def drawRect(
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
    def drawImage(
        self,
        image: Any,
        x1: Any,
        y1: Any,
        x2: Optional[Any] = ...,
        y2: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...

def getFileObject(file: Any, openFlags: str = ...): ...

class AffineMatrix:
    A: Any = ...
    def __init__(self, init: Optional[Any] = ...) -> None: ...
    def scale(self, sx: Any, sy: Any) -> None: ...
    def rotate(self, theta: Any) -> None: ...
    def translate(self, tx: Any, ty: Any) -> None: ...
