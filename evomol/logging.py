"""
Code from https://github.com/mCodingLLC/VideosSampleCode
https://github.com/mCodingLLC/VideosSampleCode/tree/master/videos/135_modern_logging
"""

import datetime as dt
import json
import logging
import logging.config
import logging.handlers
import os
import pathlib
from typing import Any

from typing_extensions import override


def init_logger(
    config_file: str = os.path.join("logging_configs", "stderr-json-file.json"),
    output_file: str = os.path.join("logs", "evomol_default.log"),
) -> logging.Logger:
    """Initialize the logger.

    Args:
        config_file (str, optional): Path to the JSON file containing the
            logging configuration settings.
            Defaults to os.path.join("logging_configs", "stderr-json-file.json").
        output_file (str, optional): Path to the output log file. Defaults to
            os.path.join("logs", "evomol_default.log").

    Returns:
        logging.Logger: The logger object.
    """

    setup_logging(
        config_file=config_file,
        output_file=output_file,
    )

    logging.basicConfig(level="INFO")
    return logging.getLogger("evomol_logger")


def get_logger() -> logging.Logger:
    """Get the logger.

    Returns:
        logging.Logger: The logger object.
    """
    return logging.getLogger("evomol_logger")


def setup_logging(
    config_file: str = os.path.join("logging_configs", "stderr-json-file.json"),
    output_file: str = "",
) -> None:
    """Setup logging configuration from a JSON file."""
    config_f = pathlib.Path(config_file)
    with open(config_f, encoding="utf-8") as f_in:
        config = json.load(f_in)

    if output_file:
        config["handlers"]["file"]["filename"] = output_file

    logging.config.dictConfig(config)


LOG_RECORD_BUILTIN_ATTRS = {
    "args",
    "asctime",
    "created",
    "exc_info",
    "exc_text",
    "filename",
    "funcName",
    "levelname",
    "levelno",
    "lineno",
    "module",
    "msecs",
    "message",
    "msg",
    "name",
    "pathname",
    "process",
    "processName",
    "relativeCreated",
    "stack_info",
    "thread",
    "threadName",
    "taskName",
}


class JSONFormatter(logging.Formatter):
    """JSON formatter for logging"""

    def __init__(
        self,
        *,
        fmt_keys: dict[str, str] | None = None,
    ):
        super().__init__()
        self.fmt_keys = fmt_keys if fmt_keys is not None else {}

    @override
    def format(self, record: logging.LogRecord) -> str:
        message = self._prepare_log_dict(record)
        return json.dumps(message, default=str)

    def _prepare_log_dict(self, record: logging.LogRecord) -> dict[str, Any]:
        always_fields = {
            "message": record.getMessage(),
            "timestamp": dt.datetime.fromtimestamp(
                record.created, tz=dt.timezone.utc
            ).isoformat(),
        }
        if record.exc_info is not None:
            always_fields["exc_info"] = self.formatException(record.exc_info)

        if record.stack_info is not None:
            always_fields["stack_info"] = self.formatStack(record.stack_info)

        message = {
            key: (
                msg_val
                if (msg_val := always_fields.pop(val, None)) is not None
                else getattr(record, val)
            )
            for key, val in self.fmt_keys.items()
        }
        message.update(always_fields)

        for key, val in record.__dict__.items():
            if key not in LOG_RECORD_BUILTIN_ATTRS:
                message[key] = val

        return message
