"""
"""

from utils import exception_handler

from typing import Tuple, List

import pysam
import os

IUPAC = {
    "A": "A",
    "T": "T",
    "C": "C",
    "G": "G",
    "a": "A",
    "t": "T",
    "c": "C",
    "g": "G",
    "N": "ATGC",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "r": "AG",
    "y": "CT",
    "s": "GC",
    "w": "AT",
    "k": "GT",
    "m": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "b": "CGT",
    "d": "AGT",
    "h": "ACT",
    "v": "ACG",
}
IUPAC_ENCODER = {
    "R": ["AG", "GA"],
    "Y": ["CT", "TC"],
    "S": ["GC", "CG"],
    "W": ["AT", "TA"],
    "K": ["GT", "TG"],
    "M": ["AC", "CA"],
    "B": ["CGT", "GCT", "TGC", "GTC", "CTG", "TCG"],
    "D": ["AGT", "GAT", "TAG", "ATG", "GTA", "TGA"],
    "H": ["ACT", "CAT", "TCA", "ATC", "CTA", "TAC"],
    "V": ["ACG", "CAG", "GAC", "AGC", "CGA", "GCA"],
    "N": ["ACGT", "CAGT", "GACT", "AGCT", "CGAT", "GCAT", "GCTA", "CGTA", "TGCA", "GTCA", "CTGA", "TCGA", "TAGC", "ATGC", "GTAC", "TGAC", "AGTC", "GATC", "CATG", "ACTG", "TCAG", "CTAG", "ATCG", "TACG"]
}
IUPAC_ENCODER_REV = {v: k for k in IUPAC_ENCODER for v in IUPAC_ENCODER[k]}

def index_fasta(fastafile: str, debug: bool) -> str:
    print(f"indexing {fastafile}")
    # fasta indexing enables fast and efficient queries 
    pysam.faidx(fastafile)
    fai = f"{fastafile}.fai"
    if not os.path.isfile(fai) or os.stat(fai).st_size <= 0:
        exception_handler(FileExistsError, f"Indexing {fastafile} failed", os.EX_SOFTWARE, debug)
    return fai

def load_genome(fastafile: str, debug: bool) -> Tuple[str, pysam.FastaFile]:
    fai = f"{fastafile}.fai"
    fasta_index = fai if os.path.isfile(fai) else index_fasta(fastafile, debug)
    try:
        fasta = pysam.FastaFile(fastafile, filepath_index=fasta_index)
        # recover sequence name (chromosome)
        assert len(fasta.references) == 1
        chrom = fasta.references[0]
    except Exception:
        exception_handler(RuntimeError, f"{fastafile} loading failed", os.EX_SOFTWARE, debug)
    return chrom, fasta


def index_vcf(vcffile: str, debug: bool) -> str:
    print(f"indexing {vcffile}")
    # vcf indexing enables faster vcf parsing and access
    pysam.tabix_index(vcffile, preset="vcf")  # index input vcf with tabix
    tbi = f"{vcffile}.tbi"
    if not os.path.isfile(tbi) or os.stat(tbi).st_size <= 0:
        exception_handler(FileExistsError, f"Indexing {vcffile} failed", os.EX_SOFTWARE, debug)
    return tbi


def load_vcf(vcffile: str, debug: bool) -> pysam.TabixFile:
    tbi = f"{vcffile}.tbi"
    tbi_index = tbi if os.path.isfile(tbi) else index_vcf(vcffile, debug)
    try:
        return pysam.TabixFile(vcffile, index=tbi_index)
    except Exception:
        exception_handler(RuntimeError, f"{vcffile} loading failed", os.EX_SOFTWARE, debug)

def encode_snp_iupac(fasta: pysam.FastaFile, chrom: str, ref: str, alt: str, pos: int, debug: bool) -> str:
    if len(ref) != 1:  # indel, likely  deletion
        return
    alleles_alt = alt.split(",")  # handle multiallelic sites
    refnt = fasta.fetch(chrom, pos - 1, pos).upper()  # 0-based position
    if refnt != ref:
       exception_handler(ValueError, f"Reference allele mismatch between FASTA and VCF file ({refnt} - {ref}, position {pos})", os.EX_SOFTWARE, debug)
    iupac_string = "".join(set([refnt] + [aa for aa in alleles_alt if len(aa) == 1]))  # recover iupac symbol
    return IUPAC_ENCODER_REV[iupac_string]

def retrieve_allele_frequency(info: List[str]):
    for field in info:
        if field.startswith("AF="):  # allele frequency
            pass

def retrieve_samples_genotype(genotype: List[str]):
    pass


def enrich_genome(fastafile: str, vcf: str, debug: bool) -> None:
    seqname, fasta = load_genome(fastafile, debug)  # load input fasta
    vcf = load_vcf(vcf, debug)  # load input vcf data
    variants_iupac = {}
    for variant in vcf.fetch():
        variant = variant.strip().split()
        chrom, snppos, snpid, vfilter, ref, alt, info = variant[:8]
        if seqname != chrom:
            exception_handler(ValueError, f"Mismatch between reference sequence name and VCF contig name ({fastafile}: {seqname} - {vcf}: {chrom})", os.EX_SOFTWARE, debug)
        if vfilter == "PASS":  # ignore variants without PASS on filter
            snppos = int(snppos)
            try:
                variants_iupac[snppos - 1]
                exception_handler(ValueError, f"Multiple variants mapped at the same position ({snppos}) in {vcf}", os.EX_DATAERR, debug)
            except KeyError:  # position not yet visited
                variants_iupac[snppos - 1] = encode_snp_iupac(fasta, chrom, ref, alt, snppos, debug)

        




if __name__ == "__main__":
    enrich_genome("test-data/chr22.fa", "test-data/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz", True)

