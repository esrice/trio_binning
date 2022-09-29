from os.path import dirname, join
from unittest.mock import patch

import pytest

from trio_binning.classify_by_kmers import main
from trio_binning.seq import readfq


def test_classify_by_kmers_help(capsys):
    with patch("sys.argv", ["classify-by-kmers", "--help"]):
        with pytest.raises(SystemExit):
            main()

    out, _ = capsys.readouterr()
    assert "Classify reads into bins" in out


def test_classify_by_kmers(capsys, tmpdir):
    with patch(
        "sys.argv",
        [
            "classify-by-kmers",
            join(dirname(__file__), "data", "test.fastq"),
            join(dirname(__file__), "data", "hapA.txt"),
            join(dirname(__file__), "data", "hapB.txt"),
            "--haplotype-a-out-prefix",
            join(tmpdir, "hapA"),
            "--haplotype-b-out-prefix",
            join(tmpdir, "hapB"),
            "--unclassified-out-prefix",
            join(tmpdir, "hapU"),
            "--no-gzip-output",
        ],
    ):
        main()

    out, _ = capsys.readouterr()
    for line in out.split("\n"):
        splits = line.strip().split("\t")
        if splits[0] == "read1":
            assert splits[1] == "U"
            assert float(splits[2]) == 0
        elif splits[0] == "readAwesome":
            assert splits[1] == "B"
            assert float(splits[3]) == 2

    hap_a_out = readfq(open(join(tmpdir, "hapA.fastq"), "r"))
    correct_out = readfq(open(join(dirname(__file__), "data", "hapA.fastq"), "r"))
    num_reads = 0
    for seq1, seq2 in zip(hap_a_out, correct_out):
        num_reads += 1
        assert seq1.name == seq2.name
        assert seq1.seq == seq2.seq
        assert seq1.qual == seq2.qual

    assert num_reads == 1
