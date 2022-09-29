"""Really really simple sequence IO

There are lots of libraries for reading and writing fastx files in
python but to my knowledge, they are all bloated and slow (e.g.,
BioPython) or read-only (e.g., screed).

This just has a single function for reading fastx files into a Read
class, which then has a print function. That's all.
"""
import gzip
import sys
from dataclasses import dataclass
from typing import Iterator, Optional, TextIO, Tuple, Union, cast


@dataclass
class Read:
    """A fastx read"""

    name: str
    """The name of the read"""
    seq: str
    """The sequence of the read"""
    qual: Optional[str] = None
    """The quality score string of the read"""

    def __str__(self):
        if self.qual:
            return f"@{self.name}\n{self.seq}\n+\n{self.qual}"
        else:
            return f">{self.name}\n{self.seq}"

    def print(self, file: TextIO = sys.stdout):
        """Print the read.

        Print the read. If it has a quality score string, print it in
        fastq format; otherwise, print it in fasta format.

        Args:
            file: the file to print the read to
        """
        print(self, file=file)


def readfq(fp: TextIO) -> Iterator[Read]:
    """Read a fastx file.

    Read a fast[aq] file, yielding a Read instance for each entry. No
    error checking is performed so it might crash. This is a lightly-
    modified version of https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for line in fp:  # search for the start of the next record
                if line[0] in ">@":  # fasta/q header line
                    last = line[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for line in fp:  # read the sequence
            if line[0] in "@+>":
                last = line[:-1]
                break
            seqs.append(line[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield Read(name, "".join(seqs), None)  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for line in fp:  # read the quality
                seqs.append(line[:-1])
                leng += len(line) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield Read(name, seq, "".join(seqs))
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield Read(name, seq, None)  # yield a fasta record instead
                break


def open_fastx_read(filename: str) -> Iterator[Read]:
    """Open a fasta/q(.gz) file for reading."""
    if filename.endswith(".gz"):
        reads = readfq(cast(TextIO, gzip.open(filename, "rt")))
    else:
        reads = readfq(open(filename, "r"))
    return reads


TextOrGzip = Union[TextIO, gzip.GzipFile]


def open_outfiles(
    haplotype_a_prefix: str,
    haplotype_b_prefix: str,
    unclassified_prefix: str,
    outfile_extension: str,
    gzip_output: bool,
) -> Tuple[TextOrGzip, TextOrGzip, TextOrGzip]:
    """Open output files based on given options.

    Args:
        haplotype_a_prefix: path prefix for haplotype A output file
        haplotype_b_prefix: path prefix for haplotype B output file
        unclassified_prefix: path prefix for unclassified output file
        outfile_extension: extension for output file (e.g., ".fa")
        gzip_output: True to gzip output files, False otherwise

    Returns:
        haplotype_a_outfile: writeable outfile for haplotype A
        haplotype_b_outfile: writeable outfile for haplotype B
        unclassified_outfile: writeable outfile for unclassified reads
    """
    haplotype_a_outfile_name = haplotype_a_prefix + outfile_extension
    haplotype_b_outfile_name = haplotype_b_prefix + outfile_extension
    unclassified_outfile_name = unclassified_prefix + outfile_extension

    haplotype_a_outfile: TextOrGzip
    haplotype_b_outfile: TextOrGzip
    unclassified_outfile: TextOrGzip

    if not gzip_output:
        haplotype_a_outfile = open(haplotype_a_outfile_name, "w")
        haplotype_b_outfile = open(haplotype_a_outfile_name, "w")
        unclassified_outfile = open(unclassified_outfile_name, "w")
    else:
        haplotype_a_outfile = gzip.open(haplotype_a_outfile_name + ".gz", "wt")
        haplotype_b_outfile = gzip.open(haplotype_b_outfile_name + ".gz", "wt")
        unclassified_outfile = gzip.open(unclassified_outfile_name + ".gz", "wt")

    return haplotype_a_outfile, haplotype_b_outfile, unclassified_outfile
