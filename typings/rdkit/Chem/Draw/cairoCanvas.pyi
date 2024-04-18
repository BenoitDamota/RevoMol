from rdkit.Chem.Draw.canvasbase import CanvasBase as CanvasBase
from typing import Any, Optional

have_cairocffi: bool
have_pango: bool
ffi: Any
defaultLibs: Any
envVar: Any
envVarSet: bool
libName: Any
libPath: Any
importError: bool
scriptPattern: Any

class Canvas(CanvasBase):
    image: Any = ...
    imageType: Any = ...
    ctx: Any = ...
    size: Any = ...
    surface: Any = ...
    fileName: Any = ...
    def __init__(
        self,
        image: Optional[Any] = ...,
        size: Optional[Any] = ...,
        ctx: Optional[Any] = ...,
        imageType: Optional[Any] = ...,
        fileName: Optional[Any] = ...,
    ) -> None: ...
    def flush(self): ...
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
    def addCircle(
        self,
        center: Any,
        radius: Any,
        color: Any = ...,
        fill: bool = ...,
        stroke: bool = ...,
        alpha: float = ...,
        **kwargs: Any,
    ) -> None: ...
