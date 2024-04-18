from qt import *
from typing import Any, Optional

class ChemdrawPanel(QWidget):
    cdx: Any = ...
    offset: int = ...
    label: Any = ...
    def __init__(
        self,
        parent: Optional[Any] = ...,
        name: str = ...,
        readOnly: int = ...,
        size: Any = ...,
    ) -> None: ...
    def pullData(self, fmt: str = ...): ...
    def setData(self, data: Any, fmt: str = ...): ...
    def resizeEvent(self, evt: Any) -> None: ...
    def __del__(self) -> None: ...
