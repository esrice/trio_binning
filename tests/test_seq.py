import os
from io import StringIO

from trio_binning import seq


def test_read_fasta():
    fasta_path = os.path.join(os.path.dirname(__file__), "data", "test.fa")

    for read in seq.readfq(open(fasta_path)):
        assert read.name.startswith("read")
        assert read.seq.startswith("G")
        assert read.seq.endswith("A")
        assert not read.qual


def test_read_fastq():
    fastq_path = os.path.join(os.path.dirname(__file__), "data", "test.fastq")

    for read in seq.readfq(open(fastq_path)):
        assert read.name.startswith("read")
        assert read.seq.startswith("G")
        assert read.seq.endswith("A")
        assert read.qual.startswith("!")


def test_write_fasta():
    sio = StringIO()
    seq.Read("read1", "AGATAGAGGACTGA").print(file=sio)
    seq.Read("read2", "AGGGGATTTTATTA").print(file=sio)
    assert sio.getvalue() == ">read1\nAGATAGAGGACTGA\n>read2\nAGGGGATTTTATTA\n"
    sio.close()


def test_write_fastq():
    sio = StringIO()
    seq.Read("read1", "AGATAGAGGACTGA", "%()%%%(%(++***").print(file=sio)
    seq.Read("read2", "AGGGGATTTTATTA", "++(*))*+%%%))(").print(file=sio)
    assert sio.getvalue() == (
        "@read1\nAGATAGAGGACTGA\n+\n%()%%%(%(++***\n"
        "@read2\nAGGGGATTTTATTA\n+\n++(*))*+%%%))(\n"
    )
    sio.close()
