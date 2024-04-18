from rdkit.rdBase import (
    AttachFileToLog as AttachFileToLog,
    DisableLog as DisableLog,
    EnableLog as EnableLog,
    LogMessage as LogMessage,
)
from typing import Any

DEBUG: int
INFO: int
WARNING: int
ERROR: int
CRITICAL: int

class logger:
    def logIt(self, dest: Any, msg: Any, *args: Any, **kwargs: Any) -> None: ...
    def debug(self, msg: Any, *args: Any, **kwargs: Any) -> None: ...
    def error(self, msg: Any, *args: Any, **kwargs: Any) -> None: ...
    def info(self, msg: Any, *args: Any, **kwargs: Any) -> None: ...
    def warning(self, msg: Any, *args: Any, **kwargs: Any) -> None: ...
    def critical(self, msg: Any, *args: Any, **kwargs: Any) -> None: ...
    def setLevel(self, val: Any) -> None: ...
