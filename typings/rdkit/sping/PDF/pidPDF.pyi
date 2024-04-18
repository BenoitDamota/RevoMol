from rdkit.sping.pid import *
from . import pdfgen as pdfgen, pdfgeom as pdfgeom, pdfmetrics as pdfmetrics
from math import ceil as ceil
from rdkit.sping import pagesizes as pagesizes
from typing import Any, Optional

DEFAULT_PAGE_SIZE: Any
font_face_map: Any
ps_font_map: Any

class PDFCanvas(Canvas):
    pdf: Any = ...
    pagesize: Any = ...
    filename: Any = ...
    drawingsize: Any = ...
    pageTransitionString: str = ...
    pageNumber: int = ...
    def __init__(
        self, size: Optional[Any] = ..., name: str = ..., pagesize: Any = ...
    ) -> None: ...
    defaultFont: Any = ...
    defaultLineColor: Any = ...
    defaultFillColor: Any = ...
    defaultLineWidth: Any = ...
    def showPage(self) -> None: ...
    def isInteractive(self): ...
    def canUpdate(self): ...
    def clear(self) -> None: ...
    def flush(self) -> None: ...
    def save(self, file: Optional[Any] = ..., format: Optional[Any] = ...) -> None: ...
    def setInfoLine(self, s: Any) -> None: ...
    def __setattr__(self, key: Any, value: Any) -> None: ...
    def resetDefaults(self) -> None: ...
    def stringWidth(self, s: Any, font: Optional[Any] = ...): ...
    def fontHeight(self, font: Optional[Any] = ...): ...
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
    def drawImage(
        self,
        image: Any,
        x1: Any,
        y1: Any,
        x2: Optional[Any] = ...,
        y2: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawLiteral(self, literal: Any) -> None: ...

def test(): ...
