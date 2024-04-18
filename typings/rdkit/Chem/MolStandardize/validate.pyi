import logging
from .errors import StopValidateError as StopValidateError
from .validations import VALIDATIONS as VALIDATIONS
from rdkit import Chem as Chem
from typing import Any

SIMPLE_FORMAT: str
LONG_FORMAT: str

class LogHandler(logging.Handler):
    logs: Any = ...
    def __init__(self) -> None: ...
    @property
    def logmessages(self): ...
    def emit(self, record: Any) -> None: ...
    def flush(self) -> None: ...
    def close(self) -> None: ...

class Validator:
    raw: Any = ...
    log: Any = ...
    handler: Any = ...
    validations: Any = ...
    def __init__(
        self,
        validations: Any = ...,
        log_format: Any = ...,
        level: Any = ...,
        stdout: bool = ...,
        raw: bool = ...,
    ) -> None: ...
    def __call__(self, mol: Any): ...
    def validate(self, mol: Any): ...

def validate_smiles(smiles: Any): ...
