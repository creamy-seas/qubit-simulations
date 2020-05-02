"""Writting output to stdout
"""

import logging


class TerminalColour:
    """Bunch of cool colours
    """

    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def setup_logging(logging_level: int):
    if logging.getLogger().handlers:
        return

    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        f"""{TerminalColour.UNDERLINE}%(funcName)s{TerminalColour.ENDC}:{TerminalColour.OKBLUE}%(lineno)d{TerminalColour.ENDC}
⮑  %(message)s"""
    )

    handler.setLevel(logging_level)
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)
