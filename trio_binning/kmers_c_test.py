from ctypes import (POINTER, Structure, byref, c_char_p, c_int, c_ubyte,
                    c_uint64, cdll)


class HashSet(Structure):
    _fields_: list = [
        ("kmers", POINTER(c_uint64)),
        ("full", POINTER(c_ubyte)),
        ("hash_size", c_int),
        ("k", c_ubyte),
        ("num_kmers", c_int),
    ]


def main():
    lib = cdll.LoadLibrary("../c/kmers.so")

    create_kmer_hash_set = lib.create_kmer_hash_set
    create_kmer_hash_set.argtypes = [c_char_p]
    create_kmer_hash_set.restype = POINTER(HashSet)

    count_kmers_in_read = lib.count_kmers_in_read
    count_kmers_in_read.argtypes = [
        c_char_p,
        c_int,
        POINTER(HashSet),
        POINTER(HashSet),
        POINTER(c_int),
        POINTER(c_int),
    ]

    hapA = create_kmer_hash_set(b"../c/hapA.txt")
    hapB = create_kmer_hash_set(b"../c/hapB.txt")

    count_a, count_b = c_int(), c_int()

    count_kmers_in_read(
        b"GAGGAGATTTAGAGTGTGAGTCGAGCATAGAGATATATA",
        hapA,
        hapB,
        byref(count_a),
        byref(count_b),
    )

    print(
        f"Found {count_a.value} kmers from hapA and {count_b.value} from hapB in read."
    )


if __name__ == "__main__":
    main()
