import xmllib
from rdkit.sping.pid import Font as Font
from typing import Any, Optional

sizedelta: int
subFraction: float
superFraction: float
greekchars: Any

class StringSegment:
    super: int = ...
    sub: int = ...
    bold: int = ...
    italic: int = ...
    underline: int = ...
    s: str = ...
    width: int = ...
    greek: int = ...
    def __init__(self) -> None: ...
    def calcNewFont(self, font: Any): ...
    def calcNewY(self, font: Any, y: Any): ...
    def dump(self) -> None: ...

class StringFormatter(xmllib.XMLParser):
    bold: int = ...
    def start_b(self, attributes: Any) -> None: ...
    def end_b(self) -> None: ...
    italic: int = ...
    def start_i(self, attributes: Any) -> None: ...
    def end_i(self) -> None: ...
    underline: int = ...
    def start_u(self, attributes: Any) -> None: ...
    def end_u(self) -> None: ...
    super: int = ...
    def start_super(self, attributes: Any) -> None: ...
    def end_super(self) -> None: ...
    sub: int = ...
    def start_sub(self, attributes: Any) -> None: ...
    def end_sub(self) -> None: ...
    greek: int = ...
    def start_greek(self, attributes: Any, letter: Any) -> None: ...
    def end_greek(self) -> None: ...
    segmentlist: Any = ...
    elements: Any = ...
    def __init__(self): ...
    def handle_data(self, data: Any) -> None: ...
    def parseSegments(self, s: Any): ...

def fontHeight(canvas: Any, font: Optional[Any] = ...): ...
def fontAscent(canvas: Any, font: Optional[Any] = ...): ...
def fontDescent(canvas: Any, font: Optional[Any] = ...): ...
def stringWidth(canvas: Any, s: Any, font: Optional[Any] = ...): ...
def rotateXY(x: Any, y: Any, theta: Any): ...
def drawString(
    canvas: Any,
    s: Any,
    x: Any,
    y: Any,
    font: Optional[Any] = ...,
    color: Optional[Any] = ...,
    angle: int = ...,
) -> None: ...
def test1() -> None: ...
def test2() -> None: ...
def allTagCombos(
    canvas: Any,
    x: Any,
    y: Any,
    font: Optional[Any] = ...,
    color: Optional[Any] = ...,
    angle: int = ...,
): ...
def stringformatTest() -> None: ...
