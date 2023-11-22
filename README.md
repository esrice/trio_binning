# Three generations haplotype binning 
This repository contains programs implementing the [trio-binning assembly method published by Koren et al.](https://www.nature.com/articles/nbt.4277)

## Installation
### Requirements
* [kmc](https://github.com/refresh-bio/KMC) (for the `find-unique-kmers` script only)
* Python
* gcc

Until I get the chance to make a conda package for this, you'll need to install
kmc seperately:
```bash
conda create -n trio_binning kmc python
conda activate trio_binning
```

Then you can download and install this package:
```bash
git clone https://github.com/esrice/trio_binning.git
cd trio_binning
pip install .
```

If you get an error like `gcc not found` you may be on a cluster where gcc
needs to be loaded explicitly. Try `module load gcc`.

## Finding unique k-mers in parental genomes
`find-unique-kmers` is a script that uses the program kmc to find k-mers that
are unique to one of the two haplotypes. Here is an example of its usage to find
unique 21-mers in short read files from a hybrid's mother and father using 8
threads:

```bash
find-unique-kmers -k 21 -p 8 mother_R1.fastq.gz,mother_R2.fastq.gz \
    father_R1.fastq.gz,father_R2.fastq.gz
```

After running this command, there will be a plain-text list of 21-mers unique to
the maternal genome in `hapA_only_kmers.txt` and the same for the paternal
genome in `hapB_only_kmers.txt`. The working directory will also contain kmc
databases for both the total set of 21-mers for each parent and the unique sets.
The numbers of unique k-mers for each haplotype is output to STDERR.

## Classifying long reads from offspring for assembly
Once you've got lists of k-mers unique to the maternal and paternal genomes,
you can use these to classify reads from the offspring into maternal and
paternal haplotypes using the program `classify-by-kmers`, like so:

```bash
classify-by-kmers \
    input_reads.fastq.gz \
    hapA_only_kmers.txt \
    hapB_only_kmers.txt \
     --haplotype-a-out-prefix classified/maternal \
     --haplotype-b-out-prefix classified/paternal \
     --unclassified-out-prefix classified/unclassified
```

This will leave you with three files in the `classified` directory:
`paternal.fastq.gz`, `maternal.fastq.gz`, and `unclassified.fastq.gz`. The
input read format is super flexible &mdash; you can give this program reads in
fasta or fastq format, gzipped or not gzipped.

## Citations
* Rice et al. (2020). "Continuous chromosome-scale haplotypes assembled from a single interspecies F1 hybrid of yak and cattle." _GigaScience_ 9(4):giaa029
* Koren et al. (2018). "Complete assembly of parental haplotypes with trio binning." _Nature Biotechnology_ 2018/10/22/online
