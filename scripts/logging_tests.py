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

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# pylint: disable=wrong-import-position, import-error

from evomol.logging import init_logger


def main() -> None:
    """Main function to test the logging configuration."""
    logger = init_logger(output_file="logs/test.log.jsonl")

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
