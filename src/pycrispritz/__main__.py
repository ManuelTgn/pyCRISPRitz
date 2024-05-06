"""
"""

from utils import sigint_handler, CRISPRITZ_COMMANDS
from pycrispritz_argparse import pyCRISPRitzArgumentParser
from addvariants import AddVariants
from version import __version__

from typing import Optional, List
from time import time, sleep

import sys

def parseargs_crispritz() -> pyCRISPRitzArgumentParser:
    parser = pyCRISPRitzArgumentParser(usage=__doc__, add_help=False)
    # general options
    group = parser.add_argument_group("Basic Options")
    group.add_argument("-h", "--help", action="help", help="Show this message and exit")
    group.add_argument("--version", action="version", help="Display pyCRISPRitz version and exit", version=__version__)
    # create  different parsers for each command
    subparsers = parser.add_subparsers(title="Commands", dest="command", description="Type 'pycrispritz <command> -h' to view the help for each command")
    # add-variants options
    parser_add_variants = subparsers.add_parser(
        f"{CRISPRITZ_COMMANDS[0]}",  # add-variants
        description="This functionality enables the integration of genetic variants "
        "into a reference genome in FASTA format, facilitating the generation "
        "of personalized or modified genomic sequences. The variants data can "
        "include single nucleotide variants (SNVs), insertions, deletions, and "
        "other genetic alterations. The variants data can include SNVs, "
        "insertions, and deletions",
        usage=f"pycrispritz {CRISPRITZ_COMMANDS[0]}",
        help="Function to incorporate variants data into a FASTA genome "
        "file",
    )
    group = parser_add_variants.add_argument_group("Options")
    group.add_argument("-g", "--genome", type=str, required=True, metavar="GENOME-DIR", help="Path to reference genome folder")
    group.add_argument("-v", "--vcf", type=str, required=True, metavar="VCF-DIR", help="Path to VCF folder")
    group.add_argument("-j", "--threads", type=int, default=1, nargs="?", metavar="NTHREADS", help="Number of threads used while adding variants. Use '0' to automatically detect and use the maximum number of available threads. Default value: %(default)s")
    return parser

def main(commandline_args: Optional[List[str]] = None) -> None:
    try:
        start = time()  # start time point
        parser = parseargs_crispritz()  # parse command line arguments
        if commandline_args is None:
            commandline_args = sys.argv[1:]
        # if no input argument print help message
        if not commandline_args:
            parser.error_noargs()
        args = parser.parse_args()  # read command line args
        assert args.command in CRISPRITZ_COMMANDS
        if args.command == CRISPRITZ_COMMANDS[0]:  # add-variants pipeline
            x = AddVariants(parser, args)
            print(x)
        print("sleeping 10 secs")
        sleep(10)
    except KeyboardInterrupt:
        sigint_handler()
    sys.stderr.write(f"\nElapsed time {(time() - start):.2f}s\n")


# ---> entry point <--- #
if __name__ == "__main__":
    main()


