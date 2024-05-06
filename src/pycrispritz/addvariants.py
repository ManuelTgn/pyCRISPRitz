"""
"""

from pycrispritz_argparse import pyCRISPRitzArgumentParser
from utils import exception_handler, CRISPRITZ_COMMANDS

from argparse import Namespace
from glob import glob

import multiprocessing
import os

class AddVariants():
    def __init__(self, parser: pyCRISPRitzArgumentParser, args: Namespace) -> None:
        self._init_genome(parser, args.genome)  # initialize genome directory
        self._init_vcf(parser, args.vcf)  # initialize vcf directory
        self._init_threads(parser, args.threads)  # initialize threads

    def __str__(self) -> str:
        assert hasattr(self, "_genome") and hasattr(self, "_vcf") and hasattr(self, "_threads")
        return (
            f"'CRISPRitz {CRISPRITZ_COMMANDS[0]}' parameters:\n"
            f"\t- Genome directory: {self._genome}\n"
            f"\t- VCF directory: {self._vcf}\n"
            f"\t- Threads number: {self._threads}\n"
        )

    def _init_genome(self, parser: pyCRISPRitzArgumentParser, genomedir: str) -> None:
        if not os.path.exists(genomedir):
            parser.error(f"Unable to locate {genomedir}")
        if not os.path.isdir(genomedir):
            parser.error(f"{genomedir} is not a directory")
        # check whether genomedir contains FASTA files
        dircontent = glob(os.path.join(genomedir, "*.fa")) + glob(os.path.join(genomedir, "*.fasta"))
        if not dircontent:
            parser.error(f"{genomedir} does not contain FASTA files ('fa' or 'fasta' extensions required)")
        self._genome = os.path.abspath(genomedir)  # avoid potential crash due to location

    def _init_vcf(self, parser: pyCRISPRitzArgumentParser, vcfdir: str) -> None:
        if not os.path.exists(vcfdir):
            parser.error(f"Unable to locate {vcfdir}")
        if not os.path.isdir(vcfdir):
            parser.error(f"{vcfdir} is not a directory")
        # check whether vcfdir contains VCF files
        dircontent = glob(os.path.join(vcfdir, "*.vcf.gz"))
        if not dircontent:
            parser.error(f"{vcfdir} does not contain VCF files ('vcf.gz' extension required)")
        self._vcf = os.path.abspath(vcfdir)  # avoid potential crash due to location
        
    def _init_threads(self, parser: pyCRISPRitzArgumentParser, threads: int) -> None:
        maxthreads = multiprocessing.cpu_count()
        if threads < 0 or threads > maxthreads:
            parser.error(f"Forbidden number of threads ({threads})")
        self._threads = maxthreads if threads == 0 else threads
