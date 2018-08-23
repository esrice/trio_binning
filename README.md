# trio_binning
This repository contains programs implementing the [trio-binning assembly method published by Koren et al.](https://www.biorxiv.org/content/early/2018/02/26/271486)

## Installation
### Requirements
* [Jellyfish](https://github.com/gmarcais/Jellyfish)
* clang (you can use `sudo apt-get install clang` to install this in ubuntu)
* liblzma (`sudo apt-get install liblzma-dev`)

## Finding unique k-mers in parental genomes
`find-unique-kmers` is a script that uses the program Jellyfish to count all the k-mers appearing in all of the paternal reads, and then automatically chooses cutoffs as described in the paper and dumps the k-mers with frequency meetings these cutoffs into output files. It is left up to the user to do the set comparison using the bash tools `sort` and `comm` as shown below, although I hope to write a program to do this automatically soon. Here is an example of how one would use these together to find all the k-mers unique to a maternal and paternal genome:

```
./find-unique-kmers -k 17 -p 10 mother_R1.fastq.gz,mother_R2.fastq.gz \
    father_R1.fastq.gz,father_R2.fastq.gz

cut -f1 haplotype0.dump | sort -T /mnt/tmp --parallel 10 > all-maternal.kmers
cut -f1 haplotype1.dump | sort -T /mnt/tmp --parallel 10 > all-paternal.kmers

comm -23 all-maternal.kmers all-paternal.kmers > maternal-only.kmers
comm -13 all-maternal.kmers all-paternal.kmers > paternal-only.kmers
```

After running these commands, there will be a plain-text list of 17-mers unique to the maternal genome in `maternal-only.kmers` and the same in for the paternal genome in `paternal-only.kmers`.

## Classifying reads from offspring
Once you've got lists of k-mers unique to the maternal and paternal genomes, you can use these to classify reads from the offspring into maternal and paternal haplotypes using the program `classify_reads`, like so:

```
./classify_reads -a maternal-only.kmers -A classified/maternal \
    -b paternal-only.kmers -B classified/paternal \
    -u offspring.fastq.gz -U classified/unclassified -c
```

This will leave you with three files in the `classified` directory: `paternal.fq.gz`, `maternal.fq.gz`, and `unclassified.fq.gz`. The input read format is super flexible &mdash; you can give this program reads in fasta or fastq format, gzipped or not gzipped, or even a bam file.

## Citation
Koren et al. (2018). "Complete assembly of parental haplotypes with trio binning." _bioRxiv_ 271486.
