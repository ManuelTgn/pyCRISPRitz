"""
A custom argument parser for pyCRISPRitz.

Handles initialization, error handling, and custom usage formatting for the tool.

Methods:
    - __init__: Initializes the pyCRISPRitzArgumentParser.
    - error: Displays an error message and exits.
    - error_noargs: Prints help message and exits.
"""


from utils import CRISPRITZ_COMMANDS
from version import __version__

from argparse import (
    SUPPRESS,
    _MutuallyExclusiveGroup,
    Action,
    ArgumentParser,
    HelpFormatter,
)
from typing import Any, Optional, Tuple, Dict, NoReturn
from colorama import Fore

import sys
import os


class pyCRISPRitzArgumentParser(ArgumentParser):
    """A custom argument parser for pyCRISPRitz.

    Args:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Attributes:
        pyCRISPRitzHelpFormatter: A custom help formatter for pyCRISPRitz.

    Methods:
        __init__: Initializes the pyCRISPRitzArgumentParser.
        error: Displays an error message and exits.
        error_noargs: Prints help message and exits.
    """

    class pyCRISPRitzHelpFormatter(HelpFormatter):
        """A custom help formatter for pyCRISPRitz.

        Args:
            usage (str): The usage string.
            actions (str): The actions string.
            groups (str): The groups string.
            prefix (Optional[str]): The prefix string. Defaults to None.

        Methods:
            add_usage: Adds usage information to the help formatter.
        """

        def add_usage(
            self, usage: str, actions: str, groups: str, prefix: Optional[str] = "None"
        ) -> None:
            """
            Add a custom usage format to the help output.

            Args:
                usage: The usage string.
                actions: Actions to include.
                groups: Groups to include.
                prefix: Optional prefix for the usage string.

            Returns:
                None
            """
            if usage != SUPPRESS:
                args = (usage, actions, groups, "")
                self._add_item(self._format_usage, args)  # define new usage format

    def __init__(self, *args: Tuple[str, Any], **kwargs: Dict[Any, Any]) -> None:
        """Initializes the pyCRISPRitzArgumentParser.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.
        """

        kwargs["formatter_class"] = self.pyCRISPRitzHelpFormatter
        # replace default usage string with tool version
        kwargs["usage"] = kwargs["usage"].replace("{version}", __version__)
        super().__init__(*args, **kwargs)

    def error(self, message: str) -> NoReturn:
        """Displays an error message and exits.

        Args:
            message (str): The error message.
        """

        # recover the command raising error
        command = sys.argv[1] if sys.argv[1] in CRISPRITZ_COMMANDS else ""
        if "invalid choice" in message:
            message = message.replace("choice", "command")  # increase clarity
        # display error message in red
        errmessage = Fore.RED + f"\nERROR: " + f"{message}\n" + Fore.RESET
        errmessage += f"\nRun pycrispritz {command} -h/--help for usage\n\n"
        sys.stderr.write(errmessage)  # write to stderr
        sys.exit(os.EX_USAGE)

    def error_noargs(self) -> NoReturn:
        """Prints help message and exits."""

        self.print_help()
        sys.exit(os.EX_USAGE)
