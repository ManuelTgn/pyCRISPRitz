"""
"""

from colorama import init, Fore
from typing import NoReturn

import sys
import os

CRISPRITZ_COMMANDS = ["add-variants"]

def sigint_handler() -> NoReturn:
    """
    Handle the SIGINT signal.

    Displays a message and exits the program with an error code.

    Returns:
        NoReturn
    """
    
    sys.stderr.write("\nCaught SIGINT. pyCRISPRitz will exit\n")
    sys.exit(os.EX_OSERR)


def exception_handler(
        exception_type: Exception, exception: str, code: int, debug: bool
) -> NoReturn:
    """
    Handle exceptions and display error messages.

    Args:
        exception_type: The type of exception.
        exception: The exception message.
        code: The error code.
        debug: A flag indicating debug mode. If set, trace full error stack.

    Returns:
        None
    """

    init()
    if debug:  # display full error stack
        raise exception_type(f"\n\n{exception}")
    # gracefully trugger runtime error and exit 
    sys.stderr.write(f"{Fore.RED}\n\nERROR: {exception}\n{Fore.RESET}")
    sys.exit(code)