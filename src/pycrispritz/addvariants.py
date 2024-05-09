"""
"""

from pycrispritz_argparse import pyCRISPRitzArgumentParser
from utils import (
    exception_handler,
    retrieve_fasta_files,
    retrieve_vcf_files,
    read_fasta_header,
    CRISPRITZ_COMMANDS,
    ENRICHED_GENOME_DIRS,
)

from typing import Dict, Tuple
from argparse import Namespace
from glob import glob

import multiprocessing
import os


class AddVariants:
    def __init__(self, parser: pyCRISPRitzArgumentParser, args: Namespace) -> None:
        self._debug = args.debug
        self._init_genome(parser, args.genome)  # initialize genome directory
        self._init_vcf(parser, args.vcf)  # initialize vcf directory
        self._init_outdir(parser, args.output)  # initialize output folder
        self._init_threads(parser, args.threads)  # initialize threads
        self._init_fasta_vcf_dict()  # initialize FASTA - VCF dictionary

    def __str__(self) -> str:
        assert (
            hasattr(self, "_genome")
            and hasattr(self, "_vcf")
            and hasattr(self, "_threads")
        )
        return (
            f"'CRISPRitz {CRISPRITZ_COMMANDS[0]}' parameters:\n"
            f"\t- Genome directory: {self._genome}\n"
            f"\t- VCF directory: {self._vcf}\n"
            f"\t- Output folder: {self._outdir}\n"
            f"\t- Threads number: {self._threads}\n"
            f"\t- Debug mode: {self._debug}\n"
        )

    def _init_genome(self, parser: pyCRISPRitzArgumentParser, genomedir: str) -> None:
        if not os.path.exists(genomedir):
            parser.error(f"Unable to locate {genomedir}")
        if not os.path.isdir(genomedir):
            parser.error(f"{genomedir} is not a directory")
        # check whether genomedir contains FASTA files
        dircontent = retrieve_fasta_files(genomedir)
        if not dircontent:
            parser.error(
                f"{genomedir} does not contain FASTA files ('fa' or 'fasta' extensions required)"
            )
        self._genome = os.path.abspath(
            genomedir
        )  # avoid potential crash due to location

    def _init_vcf(self, parser: pyCRISPRitzArgumentParser, vcfdir: str) -> None:
        if not os.path.exists(vcfdir):
            parser.error(f"Unable to locate {vcfdir}")
        if not os.path.isdir(vcfdir):
            parser.error(f"{vcfdir} is not a directory")
        # check whether vcfdir contains VCF files
        dircontent = retrieve_vcf_files(vcfdir)
        # dircontent = glob(os.path.join(vcfdir, "*.vcf.gz"))
        if not dircontent:
            parser.error(
                f"{vcfdir} does not contain VCF files ('vcf.gz' extension required)"
            )
        self._vcf = os.path.abspath(vcfdir)  # avoid potential crash due to location

    def _init_outdir(self, parser: pyCRISPRitzArgumentParser, outdir: str) -> None:
        if not os.path.exists(outdir):
            # create the output directory in current working directory
            os.mkdir(outdir)
            assert os.path.exists(outdir) and os.path.isdir(outdir)
        elif not os.path.isdir(outdir):
            parser.error(f"{outdir} is not a directory")
        # create snp and indel genome folders in outdir if not already present
        snps_genome_dir = os.path.join(outdir, ENRICHED_GENOME_DIRS[0])
        if not os.path.isdir(snps_genome_dir):
            os.mkdir(snps_genome_dir)
        indels_genome_dir = os.path.join(outdir, ENRICHED_GENOME_DIRS[1])
        if not os.path.isdir(indels_genome_dir):
            os.mkdir(indels_genome_dir)
        self._outdir = os.path.abspath(outdir)

    def _init_threads(self, parser: pyCRISPRitzArgumentParser, threads: int) -> None:
        maxthreads = multiprocessing.cpu_count()
        if threads < 0 or threads > maxthreads:
            parser.error(f"Forbidden number of threads ({threads})")
        self._threads = maxthreads if threads == 0 else threads

    def _init_fasta_vcf_dict(self) -> None:
        # recover fasta and vcf files from input directories
        fastafiles = retrieve_fasta_files(self._genome)
        vcffiles = retrieve_vcf_files(self._vcf)
        # chromosome-fasta map
        chroms_fasta_map = {read_fasta_header(f, self._debug): f for f in fastafiles}
        # chromosome-vcf map
        chroms_vcf_map = {
            chrom: next((vcf for vcf in vcffiles if chrom in vcf), None)
            for chrom in chroms_fasta_map
        }
        assert len(chroms_fasta_map) == len(chroms_vcf_map)
        # build the fasta-vcf files map
        self._fasta_vcf_map = {
            chrom: (chroms_fasta_map[chrom], chroms_vcf_map[chrom])
            for chrom in chroms_fasta_map
        }

    @property
    def fasta_vcf_dict(self) -> Dict[str, Tuple[str, str]]:
        if hasattr(self, "_fasta_vcf_map"):
            return self._fasta_vcf_map
        # always trace this error
        raise AttributeError("Fasta - VCF dictionary not yet initialized")


def enrich_genome() -> None:
    pass
