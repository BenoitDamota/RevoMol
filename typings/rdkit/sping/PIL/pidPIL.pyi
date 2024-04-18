from rdkit.sping.pid import *
from typing import Any, Optional

f: Any

class PILCanvas(Canvas):
    def __init__(self, size: Any = ..., name: str = ...) -> None: ...
    def __setattr__(self, attribute: Any, value: Any) -> None: ...
    def getImage(self): ...
    def save(self, file: Optional[Any] = ..., format: Optional[Any] = ...) -> None: ...
    def clear(self) -> None: ...
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
    def drawString(
        self,
        s: Any,
        x: Any,
        y: Any,
        font: Optional[Any] = ...,
        color: Optional[Any] = ...,
        angle: int = ...,
        **kwargs: Any,
    ): ...
    def drawImage(
        self,
        image: Any,
        x1: Any,
        y1: Any,
        x2: Optional[Any] = ...,
        y2: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...

def test(): ...
def testit(canvas: Any, s: Any, x: Any, y: Any, font: Optional[Any] = ...) -> None: ...
def test2() -> None: ...
