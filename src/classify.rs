
use kmer;
use seq::{self, fasta};
use std::result;
use std::error;
use std::fmt;

struct HaplotypeScores {
    hap_a_score: f32,
    hap_b_score: f32,
}

#[derive(Debug)]
pub enum ClassifyError {
    Kmer(kmer::KmerError),
    Reader(seq::ReaderError),
}

impl fmt::Display for ClassifyError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ClassifyError::Kmer(e) => e.fmt(f),
            ClassifyError::Reader(e) => e.fmt(f),
        }
    }
}

impl error::Error for ClassifyError {}

type Result<T> = result::Result<T, ClassifyError>;

fn count_kmers_in_read(hap_a_kmers: &kmer::KmerSet, hap_b_kmers: &kmer::KmerSet,
                       read: &seq::Record, k: usize) -> Result<(u32, u32)> {

    let mut hap_a_count: u32 = 0;
    let mut hap_b_count: u32 = 0;

    for i in 0..(read.seq().len() - k + 1) {
        let bits = kmer::get_canonical_repr(&read.seq()[i..i+k])
            .and_then(|k| kmer::kmer_to_bits(&k))
            .map_err(|e| ClassifyError::Kmer(e))?;

        if hap_a_kmers.contains(&bits) {
            hap_a_count += 1;
        }
        // hap_a_kmers and hap_b_kmers *should* be mutually exclusive sets, so
        // it shouldn't matter whether this is if or else if, but I'm doing it
        // this way so that the answer is still correct even if they aren't.
        if hap_b_kmers.contains(&bits) {
            hap_b_count += 1;
        }
    }

    Ok((hap_a_count, hap_b_count))
}


pub fn classify_unpaired(hap_a_kmers: &kmer::KmerSet,
                         hap_b_kmers: &kmer::KmerSet,
                         input_reads_filename: &str, hap_a_out_prefix: &str,
                         hap_b_out_prefix: &str, hap_u_out_prefix: &str)
    -> Result<()> {
    unimplemented!()
}

pub fn classify_paired(hap_a_kmers: &kmer::KmerSet, hap_b_kmers: &kmer::KmerSet,
                   input_reads_filename_A: &str, input_reads_filename_B: &str,
                   hap_a_out_prefix: &str, hap_b_output_prefix: &str,
                   hap_u_out_prefix: &str) -> Result<()> { unimplemented!() }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count() {
        let read = seq::Record::Fasta(fasta::Record::new(
                ">test\nACGGGCATCGCGGC")
            .unwrap());
        let hap_a_kmer_strings = ["ACGGG", "CGGGC", "AAAAA"];
        let hap_b_kmer_strings = ["ACGGG", "TTTTC", "GATAT"];
        let k: usize = 5;

        let hap_a_kmer_bits: kmer::KmerSet = hap_a_kmer_strings.iter()
            .map(|s| kmer::kmer_to_bits(s).unwrap()).collect();
        let hap_b_kmer_bits: kmer::KmerSet = hap_b_kmer_strings.iter()
            .map(|s| kmer::kmer_to_bits(s).unwrap()).collect();

        let (hap_a_counts, hap_b_counts) =
            count_kmers_in_read(&hap_a_kmer_bits, &hap_b_kmer_bits, &read, k)
            .unwrap();
        assert_eq!(hap_a_counts, 2);
        assert_eq!(hap_b_counts, 1);
    }
}
