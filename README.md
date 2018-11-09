# trio_binning
This repository contains programs implementing the [trio-binning assembly method published by Koren et al.](https://www.nature.com/articles/nbt.4277)

## Installation
### Requirements
* [kmc](https://github.com/refresh-bio/KMC)
* Python

### Binary package
A binary package compiled for 64-bit linux can be downloaded on the [release](https://github.com/esrice/trio_binning/releases) page. Just download the three executables from the latest release and you're good to go!

### From source
To compile from source, you'll need the following dependencies:
* rust (`curl https://sh.rustup.rs -sSf | sh`)
* clang (you can use `sudo apt-get install clang` to install this in ubuntu)
* liblzma (`sudo apt-get install liblzma-dev`)

Then, run the following commands to download and compile:
```
git clone https://github.com/esrice/trio_binning.git
cd trio_binning
cargo build --release
```
There will now be binaries in `target/release`.

## Finding unique k-mers in parental genomes
`find-unique-kmers` is a script that uses the program kmc to find k-mers that
are unique to one of the two haplotypes. Here is an example of its usage to find
unique 21-mers in short read files from a hybrid's mother and father using 8
threads:

```
find-unique-kmers -k 21 -p 8 mother_R1.fastq.gz,mother_R2.fastq.gz \
    father_R1.fastq.gz,father_R2.fastq.gz
```

After running this command, there will be a plain-text list of 21-mers unique to
the maternal genome in `hapA_only_kmers.txt` and the same for the paternal
genome in `hapB_only_kmers.txt`. The working directory will also contain kmc
databases for both the total set of 21-mers for each parent and the unique sets.
The numbers of unique k-mers for each haplotype is output to STDERR.

## Classifying long reads from offspring for assembly
Once you've got lists of k-mers unique to the maternal and paternal genomes, you can use these to classify reads from the offspring into maternal and paternal haplotypes using the program `classify_by_kmers`, like so:

```
classify_by_kmers -a hapA_only_kmers.txt -A classified/maternal \
    -b hapB_only_kmers.txt -B classified/paternal \
    -u offspring.fastq.gz -U classified/unclassified -c
```

This will leave you with three files in the `classified` directory: `paternal.fq.gz`, `maternal.fq.gz`, and `unclassified.fq.gz`. The input read format is super flexible &mdash; you can give this program reads in fasta or fastq format, gzipped or not gzipped, or even a bam file.

## Classifying short reads from offspring based on alignment
After you have performed the initial contig assembly, you may have short reads
that you want to bin by haplotype as well, such as Hi-C reads for scaffolding
or RNA-seq reads for annotation. Now that you have a reference genome, you no
longer need to use k-mers to determine which reads came from which haplotype;
you can instead align all reads to both haplotypes and then classify them based
on which haplotype they align to best. This package contains a program for that
as well. First, align the reads to both haplotypes using an aligner such as bwa
and sort the alignments by read name:

```
bwa mem maternal_ref.fa short_reads_R1.fq.gz short_reads_R2.fq.gz \
    | samtools view -bh - | samtools sort -n - > short_reads.maternal_ref.bam
bwa mem paternal_ref.fa short_reads_R1.fq.gz short_reads_R2.fq.gz \
    | samtools view -bh - | samtools sort -n - > short_reads.paternal_ref.bam
```

*N.B.* If you do any duplicate removal or anything else that could remove some
alignments from one bam file but not the other, do that _after_ the
classification step.

Then, use the `classify_by_alignment` program to classify each read into a
haplotype based on which assembly it aligns to best:

```
classify_by_alignment --hapA-in short_reads.maternal_ref.bam \
    --hapA-out maternal_classified.bam \
    --hapB-in short_reads.paternal_ref.bam \
    --hapB-out paternal_classified.bam
```
This command will result in a file `maternal_classified.bam` containing the
reads classified to the maternal haplotype, aligned to the maternal reference,
and a file `paternal_classified.bam` containing reads classified to the
paternal haplotype, aligned to the paternal reference. You can then post-process
these bams and use them as input to your scaffolding or annotation pipelines.

## Citation
Koren et al. (2018). "Complete assembly of parental haplotypes with trio binning." _Nature Biotechnology_ 2018/10/22/online
