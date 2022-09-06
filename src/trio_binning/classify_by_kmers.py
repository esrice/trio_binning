"""Classify reads into bins based on kmers.

This is a script for classifying sequence reads into parental bins
based on the presence of k-mers.
"""

import argparse
import gzip
from os import path
from typing import TextIO, Tuple, Union, cast

import kmers
from seq import readfq


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "reads",
        help="reads to classify into bins, in fasta/q format. Can be gzipped.",
    )
    parser.add_argument(
        "haplotype_a_kmers",
        type=kmers.create_kmer_hash_set,
        help="a list of k-mers unique to haplotype A, one per line",
    )
    parser.add_argument(
        "haplotype_b_kmers",
        type=kmers.create_kmer_hash_set,
        help="a list of k-mers unique to haplotype B, one per line",
    )
    parser.add_argument(
        "--haplotype-a-out-prefix",
        default="hapA",
        help="prefix for haplotype A output file",
    )
    parser.add_argument(
        "--haplotype-b-out-prefix",
        default="hapB",
        help="prefix for haplotype B output file",
    )
    parser.add_argument(
        "--unclassified-out-prefix",
        default="unclassified",
        help="prefix for unclassified output file",
    )
    parser.add_argument(
        "--no-gzip-output",
        action="store_true",
        help="don't gzip the output",
        default=True,
    )
    return parser.parse_args()


TextOrGzip = Union[TextIO, gzip.GzipFile]


def open_outfiles(
    haplotype_a_prefix: str,
    haplotype_b_prefix: str,
    unclassified_prefix: str,
    outfile_extension: str,
    gzip_output: bool,
) -> Tuple[TextOrGzip, TextOrGzip, TextOrGzip]:
    haplotype_a_outfile_name = haplotype_a_prefix + outfile_extension
    haplotype_b_outfile_name = haplotype_b_prefix + outfile_extension
    unclassified_outfile_name = unclassified_prefix + outfile_extension

    haplotype_a_outfile: Union[TextIO, gzip.GzipFile]
    haplotype_b_outfile: Union[TextIO, gzip.GzipFile]
    unclassified_outfile: Union[TextIO, gzip.GzipFile]

    if not gzip_output:
        haplotype_a_outfile = open(haplotype_a_outfile_name, "w")
        haplotype_b_outfile = open(haplotype_a_outfile_name, "w")
        unclassified_outfile = open(unclassified_outfile_name, "w")
    else:
        haplotype_a_outfile = gzip.open(haplotype_a_outfile_name + ".gz", "w")
        haplotype_b_outfile = gzip.open(haplotype_b_outfile_name + ".gz", "w")
        unclassified_outfile = gzip.open(unclassified_outfile_name + ".gz", "w")

    return haplotype_a_outfile, haplotype_b_outfile, unclassified_outfile


def main():
    args = parser.parse_args()

    if args.reads.endswith(".gz"):
        reads = readfq(cast(TextIO, gzip.open(args.reads, "r")))
    else:
        reads = open(args.reads, "r")

    outfile_extension = path.splitext(args.reads.rstrip(".gz"))[1]

    haplotype_a_outfile, haplotype_b_outfile, unclassified_outfile = open_outfiles(
        args.haplotype_a_out_prefix,
        args.haplotype_b_out_prefix,
        args.unclassified_out_prefix,
        outfile_extension,
        not args.no_gzip_output,
    )

    for read in reads:
        hapA_score, hapB_score = kmers.count_kmers_in_read(
            read.seq, args.haplotype_a_kmers, args.haplotype_b_kmers
        )

        # TODO be a little bit more finessed than A>B or B>A
        if hapA_score > hapB_score:
            read_bin = "A"
            read.print(file=haplotype_a_outfile)
        elif hapB_score > hapA_score:
            read_bin = "B"
            read.print(file=haplotype_b_outfile)
        else:
            read_bin = "U"
            read.print(file=unclassified_outfile)

        print("\t".join(map(str, [read.name, read_bin, hapA_score, hapB_score])))


if __name__ == "__main__":
    main()
