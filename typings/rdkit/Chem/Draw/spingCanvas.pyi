from rdkit.Chem.Draw.canvasbase import CanvasBase as CanvasBase
from rdkit.sping import pid as pid
from typing import Any, Optional

faceMap: Any

def convertColor(color: Any): ...

class Canvas(CanvasBase):
    canvas: Any = ...
    size: Any = ...
    def __init__(self, size: Any, name: Any, imageType: str = ...) -> None: ...
    def addCanvasLine(
        self,
        p1: Any,
        p2: Any,
        color: Any = ...,
        color2: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def addCanvasText(
        self, text: Any, pos: Any, font: Any, color: Any = ..., **kwargs: Any
    ): ...
    def addCanvasPolygon(
        self,
        ps: Any,
        color: Any = ...,
        fill: bool = ...,
        stroke: bool = ...,
        **kwargs: Any,
    ) -> None: ...
    def addCanvasDashedWedge(
        self,
        p1: Any,
        p2: Any,
        p3: Any,
        dash: Any = ...,
        color: Any = ...,
        color2: Optional[Any] = ...,
        **kwargs: Any,
    ) -> None: ...
    def flush(self) -> None: ...
    def save(self) -> None: ...
