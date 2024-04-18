from rdkit.sping.pid import *
from math import *
from rdkit.sping.PDF import pdfmetrics as pdfmetrics
from typing import Any, Optional

SVG_HEADER: str

class SVGCanvas(Canvas):
    size: Any = ...
    def __init__(
        self,
        size: Any = ...,
        name: str = ...,
        includeXMLHeader: bool = ...,
        extraHeaderText: str = ...,
    ) -> None: ...
    def clear(self) -> None: ...
    def flush(self) -> None: ...
    def save(self, file: Optional[Any] = ..., format: Optional[Any] = ...) -> None: ...
    def text(self): ...
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
    def drawPolygon(
        self,
        pointlist: Any,
        edgeColor: Optional[Any] = ...,
        edgeWidth: Optional[Any] = ...,
        fillColor: Any = ...,
        closed: int = ...,
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
        fillColor: Any = ...,
        dash: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def drawArc(
        self,
        x1: Any,
        y1: Any,
        x2: Any,
        y2: Any,
        theta1: int = ...,
        extent: int = ...,
        edgeColor: Optional[Any] = ...,
        edgeWidth: Optional[Any] = ...,
        fillColor: Optional[Any] = ...,
        dash: Optional[Any] = ...,
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
        fillColor: Any = ...,
        closed: int = ...,
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
    def stringWidth(self, s: Any, font: Optional[Any] = ...): ...
    def fontAscent(self, font: Optional[Any] = ...): ...
    def fontDescent(self, font: Optional[Any] = ...): ...

def test(): ...
def dashtest(): ...
def testit(canvas: Any, s: Any, x: Any, y: Any, font: Optional[Any] = ...) -> None: ...
def test2() -> None: ...
