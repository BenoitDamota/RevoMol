"""
Exemple d'utilisation de logging avec un fichier de configuration JSON.

Ce script est tiré du tuto de logging en Python :
https://github.com/mCodingLLC/VideosSampleCode/tree/master/videos/135_modern_logging

Ce code permet une sortie à la fois dans le terminal mais aussi dans un fichier
de log en JSON Lines (chaque ligne est un JSON) dans lequel on peut retrouver
toutes les informations du log et ajouter des informations supplémentaires (par
exemple, le extra={"x": "hello"}).

Il utilise aussi un système de rotation de fichier de log pour ne pas avoir un
fichier de log trop volumineux (logging.handlers.RotatingFileHandler).
"""

import json
import logging.config
import logging.handlers
import os
import pathlib
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))


def setup_logging(
    config_file: str = "logging_configs/stderr-json-file.json",
    output_file: str = "",
) -> None:
    """Setup logging configuration from a JSON file."""
    config_file = pathlib.Path(config_file)
    with open(config_file, encoding="utf-8") as f_in:
        config = json.load(f_in)

    if output_file:
        config["handlers"]["file"]["filename"] = output_file

    logging.config.dictConfig(config)


def main() -> None:
    """Main function to test the logging configuration."""
    logger = logging.getLogger("my_app")

    setup_logging(output_file="logs/test.log.jsonl")

    logging.basicConfig(level="INFO")

    logger.debug("debug message", extra={"x": "hello"})
    logger.info("info message")
    logger.warning("warning message")
    logger.error("error message")
    logger.critical("critical message")

    try:
        _ = 1 / 0
    except ZeroDivisionError:
        logger.exception("exception message")

    logger.info("coucou", extra={"x": "hello"})


if __name__ == "__main__":
    main()
