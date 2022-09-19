#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <string.h>

/*
 * Contains a hash set full of k-mers.
 */
typedef struct {
    /*
     * Contains the integer-representation of kmers
     */
    uint64_t* kmers;

    /*
     * full[i] = 1 if kmers[i] contains a kmer, 0 otherwise
     */
    unsigned char* full;

    /*
     * The size of the hash set, to be used as the length of `kmers` and `full`
     * and as the divisor of the modulo operation during lookup
     */
    int hash_size;

    /*
     * The k-mer size
     */
    unsigned char k;

    /*
     * The number of k-mers in the hash set
     */
    int num_kmers;
} hash_set;


/*
 * Convert a kmer string to a 64-bit integer representation
 *
 * Args:
 *     kmer: the kmer string to convert
 *     k: the length of the kmer (max 32)
 *
 * Returns: an integer representation of the k-mer
 */
uint64_t kmer_to_int(char* kmer, unsigned char k) {
    uint64_t kmer_int = 0;
    int i;

    for (i = 0; i < k; i++)
    {
        switch (kmer[i]) {
            case 'A':
                // no need to waste time doing 'kmer_int |= (0 << (i*2))'
                break;
            case 'C':
                kmer_int |= (1 << (i*2));
                break;
            case 'G':
                kmer_int |= (2 << (i*2));
                break;
            case 'T':
                kmer_int |= (3 << (i*2));
        }
    }

    return kmer_int;
}

void reverse_complement(char* kmer_in, char* kmer_out, unsigned char k) {
    int i;
    for (i = 0; i < k; i++)
    {
        switch (kmer_in[i]) {
            case 'A':
                kmer_out[k-i-1] = 'T';
                break;
            case 'C':
                kmer_out[k-i-1] = 'G';
                break;
            case 'G':
                kmer_out[k-i-1] = 'C';
                break;
            case 'T':
                kmer_out[k-i-1] = 'A';
                break;
        }
    }
}

/*
 * Add a k-mer to the hash.
 *
 * Args:
 *     set: the set to add the k-mer to (modifies)
 *     kmer: the k-mer string to add to the set
 */
void add_to_hash(hash_set* set, char* kmer) {
    uint64_t kmer_int = kmer_to_int(kmer, set->k);
    uint64_t position = kmer_int % set->hash_size;
    while (set->full[position])
    {
        position = (position + 1) % set->hash_size;
    }

    set->full[position] = 1;
    set->kmers[position] = kmer_int;
}

/*
 * Create a new k-mer hash set from a file full of k-mers, one per line
 */
hash_set* create_kmer_hash_set(char* kmer_file_path) {
    int num_kmers, percent_done, i;
    size_t length = 0;
    size_t characters;
    FILE* fp;

    char* line_buffer = malloc(33 * sizeof(char));
    hash_set* out_hash_set = malloc(sizeof(hash_set));

    // read through once to count number of kmers
    fp = fopen(kmer_file_path, "r");
    num_kmers = 1;
    characters = getline(&line_buffer, &length, fp);
    out_hash_set->k = (int) characters - 1;
    while (getline(&line_buffer, &length, fp) != -1) {
        num_kmers++;
    }
    out_hash_set->num_kmers = num_kmers;
    fprintf(
        stderr,
        "Found %d %d-mers in %s.\n",
        num_kmers,
        out_hash_set->k,
        kmer_file_path
    );

    // initialize hash set
    out_hash_set->hash_size = num_kmers * 4 / 3;
    //NOLINTNEXTLINE
    out_hash_set->kmers = (uint64_t*) malloc(
        out_hash_set->hash_size * sizeof(uint64_t)
    );
    out_hash_set->full = (unsigned char*) malloc(
        out_hash_set->hash_size * sizeof(unsigned char)
    );
    for (i = 0; i < out_hash_set->hash_size; i++) {
        out_hash_set->full[i] = 0;
    }

    fprintf(stderr, "Creating hash...\n");
    fp = fopen(kmer_file_path, "r");
    for (i = 0; getline(&line_buffer, &length, fp) != -1; i++) {
        add_to_hash(out_hash_set, line_buffer);
        if (i % (num_kmers/10 + 1) == 0)
        {
            percent_done = 100 * (unsigned long) i / num_kmers;
            fprintf(stderr, "%d/%d (%d%%) done\n", i, num_kmers, percent_done);
        }
    }

    fprintf(stderr, "Done!\n");
    // free up memory
    free(line_buffer);
    fclose(fp);

    return out_hash_set;
}


/*
 * Check membership of a k-mer in a hash set.
 *
 * Args:
 *     kmer: kmer string to look up
 *     empty_kmer_buffer: preallocated buffer of size sizeof(char)*(k+1) for
 *         putting reverse complements in. Just a temporary place for this
 *         function to put stuff without having to allocate the space every
 *         single time it is run, for speed reasons
 *     set: hash set in which to check for k-mer membership
 *
 * Returns: 1 if k-mer is in set, 0 otherwise
 */
char kmer_in_hash_set(char* kmer, char* empty_kmer_buffer, hash_set* set) {
    // calculate the int representations of the kmer and its reverse
    // complement, looking up only the lesser of the two
    uint64_t kmer_int = kmer_to_int(kmer, set->k);
    reverse_complement(kmer, empty_kmer_buffer, set->k);
    uint64_t kmer_revcomp_int = kmer_to_int(empty_kmer_buffer, set->k);

    kmer_int = kmer_int < kmer_revcomp_int ? kmer_int : kmer_revcomp_int;

    int position = kmer_int % set->hash_size;
    while (set->full[position])
    {
        if (set->kmers[position] == kmer_int)
        {
            return 1;
        }
        position = (position + 1) % set->hash_size;
    }

    return 0;
}

void count_kmers_in_read(
    char* read,
    hash_set* haplotype_A,
    hash_set* haplotype_B,
    int* count_A,
    int* count_B
) {
    int i, j;
    char* kmer = malloc((haplotype_A->k + 1) * sizeof(char));
    char* kmer_revcomp = malloc((haplotype_A->k + 1) * sizeof(char));
    kmer[haplotype_A->k] = '\0';

    *count_A = 0;
    *count_B = 0;

    int read_length = strlen(read);

    for (i = 0; i < read_length - haplotype_A->k + 1; i++)
    {
        for (j = 0; j < haplotype_A->k; j++)
            kmer[j] = read[i+j];
        if (kmer_in_hash_set(kmer, kmer_revcomp, haplotype_A))
            (*count_A)++;
        if (kmer_in_hash_set(kmer, kmer_revcomp, haplotype_B))
            (*count_B)++;
    }

    free(kmer);
    free(kmer_revcomp);
}

int main() {
    fprintf(stderr, "Reading hapA k-mers into hash...\n");
    hash_set* hapA = create_kmer_hash_set("hapA.txt");
    fprintf(
        stderr,
        "Finished reading %d-mers into hash of size %d.\n",
        hapA->k,
        hapA->hash_size
    );

    fprintf(stderr, "Reading hapB k-mers into hash...\n");
    hash_set* hapB = create_kmer_hash_set("hapB.txt");
    fprintf(
        stderr,
        "Finished reading %d-mers into hash of size %d.\n",
        hapA->k,
        hapA->hash_size
    );

    int count_A = 0, count_B = 0;
    count_kmers_in_read(
        "GAGGAGATTTAGAGTGTGAGTCGAGCATAGAGATATATA",
        hapA,
        hapB,
        &count_A,
        &count_B
    );
    fprintf(
        stderr,
        "Found %d k-mers from hapA and %d from hapB.\n",
        count_A,
        count_B
    );

    free(hapA);
    free(hapB);
    return 0;
}
