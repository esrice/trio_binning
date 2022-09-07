import os.path

import kmers


def test_build_kmer_hash():
    kmers_path = os.path.join(os.path.dirname(__file__), "data", "hapA.txt")
    kmer_set = kmers.create_kmer_hash_set(kmers_path)
    assert kmers.get_number_kmers_in_set(kmer_set) == 6


def test_count_kmers_in_read():
    hap_a_path = os.path.join(os.path.dirname(__file__), "data", "hapA.txt")
    hap_b_path = os.path.join(os.path.dirname(__file__), "data", "hapB.txt")

    hap_a_set = kmers.create_kmer_hash_set(hap_a_path)
    hap_b_set = kmers.create_kmer_hash_set(hap_b_path)

    assert kmers.count_kmers_in_read(
        "GAGGAGATTTAGAGTGTGAGTCGAGCATAGAGATATATA", hap_a_set, hap_b_set
    ) == (1, 2)
