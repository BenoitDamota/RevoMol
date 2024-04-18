from typing import Any

class MolSupplier:
    def __init__(self) -> None: ...
    def Reset(self) -> None: ...
    def __iter__(self) -> Any: ...
    def next(self): ...
    def NextMol(self) -> None: ...
    __next__: Any = ...
