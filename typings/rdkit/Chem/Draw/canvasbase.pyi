from typing import Any, Optional

class CanvasBase:
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
    ) -> None: ...
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
    def flush(self) -> None: ...
