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


class QuantumLogger:
    """Pretty print logging.debug(), logging.warn() and warning.err()
    """

    def __init__(self, logging_level):

        # Remove any loggers already present
        if len(logging.getLogger().handlers) > 0:
            logging.getLogger().handlers = []

        formatter = logging.Formatter(
            f"""{TerminalColour.UNDERLINE}%(module)s:{TerminalColour.OKBLUE}%(lineno)d{TerminalColour.ENDC}:{TerminalColour.UNDERLINE}%(funcName)s{TerminalColour.ENDC}
⮑  %(message)s
"""
        )

        handler = logging.StreamHandler()
        handler.setLevel(logging_level)
        handler.setFormatter(formatter)

        logging.getLogger().setLevel(logging.DEBUG)
        logging.getLogger().addHandler(handler)
