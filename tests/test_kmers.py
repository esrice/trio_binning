import os.path

import pytest

from trio_binning import kmers


def test_build_kmer_hash():
    kmers_path = os.path.join(os.path.dirname(__file__), "data", "hapA.txt")
    kmer_set = kmers.create_kmer_hash_set(kmers_path)
    assert kmers.get_number_kmers_in_set(kmer_set) == 4


def test_count_kmers_in_read():
    hap_a_path = os.path.join(os.path.dirname(__file__), "data", "hapA.txt")
    hap_b_path = os.path.join(os.path.dirname(__file__), "data", "hapB.txt")

    hap_a_set = kmers.create_kmer_hash_set(hap_a_path)
    hap_b_set = kmers.create_kmer_hash_set(hap_b_path)

    assert kmers.count_kmers_in_read(
        "CTTATCATGTCTTTGTTTTCAAAGCTTCTTAGAGGTTTTTTTTTTTGGTGTTAATTGGCATAAATTATGGCT",
        hap_a_set,
        hap_b_set,
    ) == (2, 1)


@pytest.mark.parametrize(
    "kmer_str,kmer_int",
    [
        ("ATGCTAGCTAGAGAGAGAGGA", 696357446508),
        ("TTTTTTTTTTTTTTTTTTTTTTTTTTTT", 72057594037927935),
        ("TTTTTTTTTTTTTTTAGGCCCACTTTTT", 72006304612220927),
        ("AAAAAAAAAAAAAAAAAAAAAAAAAA", 0),
        ("GGGAGGGAGGGAGGGAGGGAGGGAGGG", 11868309606246954),
    ],
)
def test_kmer_to_int(kmer_str, kmer_int):
    assert kmers.kmer_to_int(kmer_str) == kmer_int


@pytest.mark.parametrize(
    "kmer,revcomp_kmer",
    [
        ("ATGCTAGCTAGAGAGAGAGGA", "TCCTCTCTCTCTAGCTAGCAT"),
        ("TTTTTTTTTTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        ("TTTTTTTTTTTTTTTAGGCCCACTTTTT", "AAAAAGTGGGCCTAAAAAAAAAAAAAAA"),
        ("AAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTT"),
        ("GGGAGGGAGGGAGGGAGGGAGGGAGGG", "CCCTCCCTCCCTCCCTCCCTCCCTCCC"),
    ],
)
def test_reverse_complement(kmer, revcomp_kmer):
    assert kmers.reverse_complement(kmer) == revcomp_kmer
