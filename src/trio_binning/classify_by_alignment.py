"""Bin reads based on alignment.

Classify reads from an F1 hybrid into parental bins based on whether
each one aligns better to the maternal or paternal assembly.
"""
import argparse
from typing import Iterator

import mappy
import pysam

from trio_binning import seq


def open_sam_or_bam(filename: str) -> pysam.AlignmentFile:
    """Open a sam or bam file with correct mode

    Open a sam or bam file. If it is a sam file, open it in "r" mode.
    If it is a bam file, open it in "b" mode.

    Args:
        filename: path to the sam or bam file

    Returns: the opened file
    """
    open_mode = "r"
    if filename.endswith(".bam"):
        open_mode += "b"
    return pysam.AlignmentFile(filename, open_mode)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "haplotype_a_assembly",
        help="assembly of haplotype A for alignment",
    )
    parser.add_argument(
        "haplotype_b_assembly",
        help="assembly of haplotype B for alignment",
    )
    parser.add_argument(
        "reads_to_classify",
        type=seq.open_fastx_read,
        help="reads to classify into bins, in fasta/q format",
    )
    parser.add_argument(
        "-m",
        "--minimap-preset",
        default="sr",
        help="minimap preset mode (see minimap2 documentation for options)",
    )
    parser.add_argument(
        "-2",
        "--reads-to-classify-r2",
        type=seq.open_fastx_read,
        help="R2 reads to classify",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    aligner_hap_a = mappy.Aligner(args.haplotype_a_assembly, preset=args.minimap_preset)
    aligner_hap_b = mappy.Aligner(args.haplotype_b_assembly, preset=args.minimap_preset)

    hap_a_out, hap_b_out, hap_u_out = seq.open_outfiles()

    if args.reads_to_classify_r2:
        classify_paired(
            aligner_hap_a,
            aligner_hap_b,
            args.reads_to_classify,
            args.reads_to_classify_r2,
        )
    else:
        classify_single(
            aligner_hap_a,
            aligner_hap_b,
            args.reads_to_classify,
            hap_a_out,
            hap_b_out,
            hap_u_out,
        )


def classify_single(
    aligner_hap_a: mappy.Aligner,
    aligner_hap_b: mappy.Aligner,
    reads_to_classify: Iterator[seq.Read],
    hap_a_out: seq.TextOrGzip,
    hap_b_out: seq.TextOrGzip,
    hap_u_out: seq.TextOrGzip,
):
    for read in reads_to_classify:
        alignment_a = max(aligner_hap_a.align(read.seq), key=lambda a: a.mlen)
        alignment_b = max(aligner_hap_b.align(read.seq), key=lambda a: a.mlen)
        if alignment_a.mlen > alignment_b.mlen:
            print(read, file=hap_a_out)
        elif alignment_a.mlen < alignment_b.mlen:
            print(read, file=hap_b_out)
        else:
            print(read, file=hap_u_out)


def classify_paired(
    aligner_hap_a, aligner_hap_b, reads_to_classify_r1, reads_to_classify_r2
):
    pass


if __name__ == "__main__":
    main()
