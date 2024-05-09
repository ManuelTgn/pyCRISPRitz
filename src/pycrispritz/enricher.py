"""
"""

from utils import exception_handler, ENRICHED_GENOME_DIRS

from typing import Tuple, List, Dict

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
    "N": [
        "ACGT",
        "CAGT",
        "GACT",
        "AGCT",
        "CGAT",
        "GCAT",
        "GCTA",
        "CGTA",
        "TGCA",
        "GTCA",
        "CTGA",
        "TCGA",
        "TAGC",
        "ATGC",
        "GTAC",
        "TGAC",
        "AGTC",
        "GATC",
        "CATG",
        "ACTG",
        "TCAG",
        "CTAG",
        "ATCG",
        "TACG",
    ],
}
IUPAC_ENCODER_REV = {v: k for k in IUPAC_ENCODER for v in IUPAC_ENCODER[k]}


def index_fasta(fastafile: str, debug: bool) -> str:
    """
    Index a FASTA file with pysam for efficient queries.

    Args:
        fastafile (str): The path to the FASTA file to index.
        debug (bool): A flag indicating debug mode.

    Returns:
        str: The path to the generated index file.
    Raises:
        FileExistsError: If indexing the FASTA file fails.
    """

    # fasta indexing enables fast and efficient queries
    pysam.faidx(fastafile)
    fai = f"{fastafile}.fai"
    if not os.path.isfile(fai) or os.stat(fai).st_size <= 0:
        exception_handler(
            FileExistsError, f"Indexing {fastafile} failed", os.EX_SOFTWARE, debug
        )
    return fai


def load_fasta(fastafile: str, debug: bool) -> Tuple[str, pysam.FastaFile]:
    """
    Load a FASTA file and its corresponding index.

    Args:
        fastafile (str): The path to the FASTA file to load.
        debug (bool): A flag indicating debug mode.

    Returns:
        Tuple[str, pysam.FastaFile]: A tuple containing the chromosome name and
            the loaded FASTA file.
    Raises:
        RuntimeError: If loading the FASTA file fails.
    """

    fai = f"{fastafile}.fai"
    fasta_index = fai if os.path.isfile(fai) else index_fasta(fastafile, debug)
    try:
        fasta = pysam.FastaFile(fastafile, filepath_index=fasta_index)
        # recover sequence name (chromosome)
        assert len(fasta.references) == 1
        chrom = fasta.references[0]
    except Exception:
        exception_handler(
            RuntimeError, f"{fastafile} loading failed", os.EX_SOFTWARE, debug
        )
    return chrom, fasta


def index_vcf(vcffile: str, debug: bool) -> str:
    """
    Index a VCF file with tabix for faster parsing and access.

    Args:
        vcffile (str): The path to the VCF file to index.
        debug (bool): A flag indicating debug mode.

    Returns:
        str: The path to the generated index file.
    Raises:
        FileExistsError: If indexing the VCF file fails.
    """

    print(f"indexing {vcffile}")
    # vcf indexing enables faster vcf parsing and access
    pysam.tabix_index(vcffile, preset="vcf")  # index input vcf with tabix
    tbi = f"{vcffile}.tbi"
    if not os.path.isfile(tbi) or os.stat(tbi).st_size <= 0:
        exception_handler(
            FileExistsError, f"Indexing {vcffile} failed", os.EX_SOFTWARE, debug
        )
    return tbi


def load_vcf(vcffile: str, debug: bool) -> Tuple[pysam.TabixFile, List[str]]:
    """
    Load a VCF file and its corresponding index for efficient access.

    Args:
        vcffile (str): The path to the VCF file to load.
        debug (bool): A flag indicating debug mode.

    Returns:
        Tuple[pysam.TabixFile, List[str]]: A tuple containing the loaded VCF
            file and a list of samples extracted from the VCF file.
    Raises:
        RuntimeError: If loading the VCF file fails.
    """

    tbi = f"{vcffile}.tbi"
    tbi_index = tbi if os.path.isfile(tbi) else index_vcf(vcffile, debug)
    try:
        vcf = pysam.TabixFile(vcffile, index=tbi_index)
        # recover vcf samples (from 9th col)
        samples = vcf.header[-1].strip().split()[9:]
        return vcf, samples
    except Exception:
        exception_handler(
            RuntimeError, f"{vcffile} loading failed", os.EX_SOFTWARE, debug
        )


def encode_snp_iupac(
    fasta: pysam.FastaFile, chrom: str, ref: str, alt: str, pos: int, debug: bool
) -> str:
    if len(ref) != 1:  # indel, likely  deletion
        return
    alleles_alt = alt.split(",")  # handle multiallelic sites
    refnt = fasta.fetch(chrom, pos - 1, pos).upper()  # 0-based position
    if refnt != ref:
        exception_handler(
            ValueError,
            f"Reference allele mismatch between FASTA and VCF file ({refnt} - {ref}, position {pos})",
            os.EX_SOFTWARE,
            debug,
        )
    iupac_string = "".join(
        set([refnt] + [aa for aa in alleles_alt if len(aa) == 1])
    )  # recover iupac symbol
    return IUPAC_ENCODER_REV[iupac_string]


def retrieve_allele_frequency(info: List[str]):
    for field in info:
        if field.startswith("AF="):  # allele frequency
            pass


def retrieve_samples_genotype(
    genotypes: List[str], samples: List[str], vcffile: str, debug: bool
):
    try:
        return [f"{samples[i]}:{gt}" for i, gt in enumerate(genotypes) if "1" in gt]
    except Exception:
        exception_handler(
            RuntimeError,
            f"An error occurred while retrieving samples genotypes in {vcffile}",
            os.EX_SOFTWARE,
            debug,
        )


def build_snps_dict(
    chrom: str, snppos: int, genotypes: List[str], samples: List[str], debug: bool
):
    samples_genotypes = retrieve_samples_genotype(genotypes, samples)


def build_iupac_variants_dict(
    fasta: pysam.FastaFile,
    fastafile: str,
    vcf: pysam.TabixFile,
    vcffile: str,
    seqname: str,
    debug: bool,
) -> Dict[int, str]:
    variants_iupac = {}
    for variant in vcf.fetch():
        variant = variant.strip().split()
        chrom, snppos, snpid, ref, alt, _, vfilter, info = variant[:8]
        genotypes = variant[9:]
        if seqname != chrom:
            exception_handler(
                ValueError,
                f"Mismatch between reference sequence name and VCF contig name ({fastafile}: {seqname} - {vcffile}: {chrom})",
                os.EX_SOFTWARE,
                debug,
            )
        if vfilter == "PASS":  # ignore variants without PASS on filter
            snppos = int(snppos)
            try:
                _ = variants_iupac[snppos - 1]
                exception_handler(
                    ValueError,
                    f"Multiple variants mapped at the same position ({snppos}) in {vcf}",
                    os.EX_DATAERR,
                    debug,
                )
            except KeyError:  # position not yet visited
                variants_iupac[snppos - 1] = encode_snp_iupac(
                    fasta, chrom, ref, alt, snppos, debug
                )
        break
    return variants_iupac


def enrich_sequence(
    fasta: pysam.FastaFile,
    variants_iupac: Dict[int, str],
    chrom: str,
    outdir: str,
    debug: bool,
):
    sequence_enr = ""
    start = 0
    for pos in variants_iupac:
        sequence_enr += fasta.fetch(chrom, start, pos).upper() + variants_iupac[pos]
        start = pos + 1
    sequence_enr += fasta.fetch(chrom, start, fasta.lengths[0]).upper()
    fastaenriched, ext = os.path.splitext(os.path.basename(fasta.filename.decode()))
    fastaenriched = os.path.join(
        outdir, ENRICHED_GENOME_DIRS[0], f"{fastaenriched}.enriched{ext}"
    )
    try:
        with open(fastaenriched, mode="w") as outfile:
            outfile.write(f">{chrom}\n{sequence_enr}\n")
    except IOError:
        exception_handler(
            IOError,
            f"An error occurred while writing enriched sequence {fastaenriched}",
            os.EX_CANTCREAT,
            debug,
        )

    with open(fastaenriched) as infile:
        infile.readline()
        seqnew = infile.read()

    with open("SNPs_genome/test-output_enriched/chr22.enriched.fa") as infile:
        infile.readline()
        seqold = infile.read()

    print(len(seqnew), len(seqold))
    print(seqnew == seqold)
    print(seqnew[pos], seqold[pos])


def enrich(fastafile: str, vcffile: str, outdir: str, debug: bool) -> None:
    seqname, fasta = load_fasta(fastafile, debug)  # load input fasta
    vcf, samples = load_vcf(vcffile, debug)  # load input vcf data
    variants_iupac = build_iupac_variants_dict(
        fasta, fastafile, vcf, vcffile, seqname, debug
    )
    enrich_sequence(fasta, variants_iupac, seqname, outdir, debug)


if __name__ == "__main__":
    enrich(
        "test-data/chr22.fa",
        "test-data/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz",
        "test-output",
        True,
    )
