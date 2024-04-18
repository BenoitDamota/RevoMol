from rdkit.Chem.Draw.canvasbase import CanvasBase as CanvasBase
from typing import Any, Optional

class Canvas(CanvasBase):
    size: Any = ...
    def __init__(self, size: Any, name: str = ..., imageType: str = ...) -> None: ...
    def rescalePt(self, p1: Any): ...
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
    def addCanvasPolygon(self, ps: Any, color: Any = ..., **kwargs: Any) -> None: ...
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
