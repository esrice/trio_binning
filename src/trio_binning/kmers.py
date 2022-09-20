"""Functions for counting k-mers quickly.

This module contains functions that interface with the C library for
counting k-mers, in order to hide all the ugly C pointer stuff from
python code that just wants to count k-mers.

>>> hapA = kmers.create_kmer_hash_set("../c/hapA.txt")
Found 7 13-mers in file.
>>> hapB = kmers.create_kmer_hash_set("../c/hapB.txt")
Found 6 13-mers in file.
>>> kmers.count_kmers_in_read("GAGGAGATTTAGAGTGTGAGTCGAGCATAGAGATATATA", hapA, hapB)
(1, 2)
"""
import sys
from ctypes import (
    POINTER,
    Structure,
    byref,
    c_char_p,
    c_int,
    c_ubyte,
    c_uint64,
    cdll,
    pointer,
)
from importlib.machinery import EXTENSION_SUFFIXES
from os.path import dirname, isfile, join
from typing import TYPE_CHECKING, Tuple

correct_library_file = ""
for extension in EXTENSION_SUFFIXES:
    possible_library_file = join(dirname(__file__), "kmers_c" + extension)
    if isfile(possible_library_file):
        correct_library_file = possible_library_file

if correct_library_file == "":
    raise ImportError("Cannot load kmers_c library. Is it installed correctly?")
lib = cdll.LoadLibrary(correct_library_file)


class _HashSet(Structure):
    """Container for a c struct containing a k-mer hash set

    This is just a container for a c struct. It is necessary to define
    it here in order to allow python code interfacing with code taking
    a pointer to this struct as an object, or returning a pointer to it
    as an argument, to understand what is getting passed or returned.

    `HashSet` (no underscore) is a pointer to this class, and the
    actual type that gets accepted and returned by the C functions.
    """

    _fields_: list = [
        ("kmers", POINTER(c_uint64)),
        ("full", POINTER(c_ubyte)),
        ("hash_size", c_int),
        ("k", c_ubyte),
        ("num_kmers", c_int),
    ]


create_kmer_hash_set_c = lib.create_kmer_hash_set
create_kmer_hash_set_c.argtypes = [c_char_p]
create_kmer_hash_set_c.restype = POINTER(_HashSet)

count_kmers_in_read_c = lib.count_kmers_in_read
count_kmers_in_read_c.argtypes = [
    c_char_p,
    POINTER(_HashSet),
    POINTER(_HashSet),
    POINTER(c_int),
    POINTER(c_int),
]

# this is ugly as sin, but necessary because mypy is ok with `pointer`
# as a subscriptable type while python runtime is not
if TYPE_CHECKING:
    HashSet = pointer[_HashSet]
else:
    HashSet = pointer


def create_kmer_hash_set(kmer_file_path: str) -> HashSet:
    """Read a list of k-mers into a set.

    Reads a list of k-mers into a quickly searchable hash set.

    Args:
        kmer_file_path: the path to a file containing a set of k-mers,
            one per line.

    Returns:
        a quickly searchable set of these k-mers that can be passed to
        `count_kmers_in_read`
    """
    if not isfile(kmer_file_path):
        raise IOError(f"Specified file {kmer_file_path} does not exist or is not file.")

    print(f"Reading k-mers in {kmer_file_path}...", file=sys.stderr)

    return create_kmer_hash_set_c(kmer_file_path.encode("utf-8"))


def count_kmers_in_read(
    read: str, kmers_hap_a: HashSet, kmers_hap_b: HashSet
) -> Tuple[int, int]:
    """Count k-mers in read and two sets

    Counts the k-mers in a read, keeping track of how many are in two
    different sets, respectively.

    Args:
        read: a string containing a DNA sequence read. May only contain
            [ACGT].
        kmers_hap_a: a hash set containing all k-mers in haplotype A
        kmers_hap_b: a hash set containing all k-mers in haplotype B

    Returns:
        A tuple of two ints, where the first is the number of k-mers in
        the read found in the haplotype A set, and the second is the
        number of k-mers in the read found in the haplotype B set
    """
    count_a, count_b = c_int(), c_int()

    count_kmers_in_read_c(
        read.encode("utf-8"),
        kmers_hap_a,
        kmers_hap_b,
        byref(count_a),
        byref(count_b),
    )

    return count_a.value, count_b.value


def get_number_kmers_in_set(kmer_hash_set: HashSet) -> int:
    """Look up the number of k-mers in a hash set"""
    return kmer_hash_set.contents.num_kmers
