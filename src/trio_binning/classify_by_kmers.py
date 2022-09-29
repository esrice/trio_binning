"""Classify reads into bins based on kmers.

This is a script for classifying sequence reads into parental bins
based on the presence of k-mers.
"""

import argparse
from os import path
from typing import Tuple

from trio_binning import kmers, seq


def parse_args():
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
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
        default=False,
    )
    return parser.parse_args()


def calculate_scaling_factors(
    haplotype_a_kmers: kmers.HashSet, haplotype_b_kmers: kmers.HashSet
) -> Tuple[float, float]:
    """Calculate the scaling factors for k-mer scores

    Args:
        haplotype_a_kmers: set of k-mers in haplotype A
        haplotype_b_kmers: set of k-mers in haplotype B

    Returns:
        scaling_factor_a: scaling factor by which haplotype A counts
            should be multiplied
        scaling_factor_b: scaling factor by which haplotype B counts
            should be multiplied
    """
    num_kmers_a = kmers.get_number_kmers_in_set(haplotype_a_kmers)
    num_kmers_b = kmers.get_number_kmers_in_set(haplotype_b_kmers)
    max_num_kmers = max(num_kmers_a, num_kmers_b)
    scaling_factor_a = 1.0 * max_num_kmers / num_kmers_a
    scaling_factor_b = 1.0 * max_num_kmers / num_kmers_b
    return scaling_factor_a, scaling_factor_b


def main():
    """Main method of program"""
    args = parse_args()

    reads = seq.open_fastx_read(args.reads)

    haplotype_a_outfile, haplotype_b_outfile, unclassified_outfile = seq.open_outfiles(
        args.haplotype_a_out_prefix,
        args.haplotype_b_out_prefix,
        args.unclassified_out_prefix,
        path.splitext(args.reads.rstrip(".gz"))[1],
        not args.no_gzip_output,
    )

    scaling_factor_a, scaling_factor_b = calculate_scaling_factors(
        args.haplotype_a_kmers,
        args.haplotype_b_kmers,
    )

    for read in reads:
        hap_a_count, hap_b_count = kmers.count_kmers_in_read(
            read.seq, args.haplotype_a_kmers, args.haplotype_b_kmers
        )

        hap_a_score = hap_a_count * scaling_factor_a
        hap_b_score = hap_b_count * scaling_factor_b

        if hap_a_score > hap_b_score:
            read_bin = "A"
            read.print(file=haplotype_a_outfile)
        elif hap_b_score > hap_a_score:
            read_bin = "B"
            read.print(file=haplotype_b_outfile)
        else:
            read_bin = "U"
            read.print(file=unclassified_outfile)

        print("\t".join(map(str, [read.name, read_bin, hap_a_score, hap_b_score])))


if __name__ == "__main__":
    main()
