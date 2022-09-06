"""Really really simple sequence IO

There are lots of libraries for reading and writing fastx files in
python but to my knowledge, they are all bloated and slow (e.g.,
BioPython) or read-only (e.g., screed).

This just has a single function for reading fastx files into a Read
class, which then has a print function. That's all.
"""
import sys
from abc import ABC
from dataclasses import dataclass
from typing import IO, Iterator, Optional, TextIO


@dataclass
class Read:
    """A fastx read"""

    name: str
    """The name of the read"""
    seq: str
    """The sequence of the read"""
    qual: Optional[str]
    """The quality score string of the read"""

    def print(self, file: TextIO = sys.stdout):
        """Print the read.

        Print the read. If it has a quality score string, print it in
        fastq format; otherwise, print it in fasta format.

        Args:
            file: the file to print the read to
        """
        if self.qual:
            print(f"@{self.name}\n{self.seq}\n+\n{self.qual}", file=file)
        else:
            print(f">{self.name}\n{self.seq}", file=file)


def readfq(fp: TextIO) -> Iterator[Read]:
    """Read a fastx file.

    Read a fast[aq] file, yielding a Read instance for each entry. No
    error checking is performed so it might crash. This is a lightly-
    modified version of https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield Read(name, "".join(seqs), None)  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield Read(name, seq, "".join(seqs))
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield Read(name, seq, None)  # yield a fasta record instead
                break
