"""
"""

from colorama import init, Fore
from typing import NoReturn, List
from glob import glob

import sys
import os

CRISPRITZ_COMMANDS = ["add-variants"]
FASTAEXT = ["fa", "fasta"]
ENRICHED_GENOME_DIRS = ["SNPs_genome", "INDELs_genome"]

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

def retrieve_fasta_files(path: str) -> List[str]:
    """
    Retrieve a list of paths to FASTA files in the specified directory.

    Args:
        path: The directory path to search for FASTA files.

    Returns:
        List[str]: A list of paths to FASTA files found in the specified directory.
    """

    return [
        f for ext in FASTAEXT for f in glob(os.path.join(path, f"*.{ext}"))
    ]

def retrieve_vcf_files(path: str) -> List[str]:
    """
    Retrieve a list of paths to VCF files in the specified directory.

    Args:
        path: The directory path to search for VCF files.

    Returns:
        List[str]: A list of paths to VCF files found in the specified directory.
    """

    return list(glob(os.path.join(path, "*.vcf.gz")))

def read_fasta_header(fastafile: str, debug: bool) -> str:
    try:
        with open(fastafile, mode="r") as infile:
            for line in infile:
                if line.startswith(">"):
                    return line[1:].strip()
    except IOError:
        exception_handler(IOError, f"An error occurred while recovering chromosome name in {fastafile}", os.EX_OSFILE, debug)