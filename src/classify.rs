use kmer;
use seq;
use std::result;
use std::error;
use std::fmt;
use std::cmp;

struct HaplotypeScores {
    hap_a_score: f32,
    hap_b_score: f32,
}

type Result<T> = result::Result<T, Box<dyn error::Error>>;

fn count_kmers_in_read(hap_a_kmers: &kmer::KmerSet, hap_b_kmers: &kmer::KmerSet,
                       read: &seq::SeqRecord, k: usize) -> Result<(u32, u32)> {

    let mut hap_a_count: u32 = 0;
    let mut hap_b_count: u32 = 0;

    for i in 0..(read.seq.len() - k + 1) {
        let bits = kmer::get_canonical_repr(&read.seq[i..i+k])
            .and_then(|k| kmer::kmer_to_bits(&k))?;

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

/// look at the sizes of the k-mer sets and use these to calculating scaling
/// factors. The original program divides read haplotype counts by the size
/// of the k-mer set for that haplotype, resulting in a really tiny number.
/// In order to make it more pleasant for humans to look at, we multiply
/// both scaling factors by the size of the larger k-mer set so that the
/// scaling factors are close to 1.
pub fn calc_scaling_factors(hap_a_kmers: &kmer::KmerSet,
                            hap_b_kmers: &kmer::KmerSet) -> (f32, f32) {

    let num_hap_a_kmers = hap_a_kmers.len();
    let num_hap_b_kmers = hap_b_kmers.len();
    let max_num_kmers = cmp::max(num_hap_a_kmers, num_hap_b_kmers);
    let scaling_factor_a = (max_num_kmers as f32) / (num_hap_a_kmers as f32);
    let scaling_factor_b = (max_num_kmers as f32) / (num_hap_b_kmers as f32);

    (scaling_factor_a, scaling_factor_b)
}

pub fn classify_unpaired(hap_a_kmers: &kmer::KmerSet,
                         hap_b_kmers: &kmer::KmerSet,
                         input_reads_filename: &str, hap_a_out_prefix: &str,
                         hap_b_out_prefix: &str, hap_u_out_prefix: &str)
    -> Result<()> {

    // calculate read-count scaling factors
    let (scaling_factor_a, scaling_factor_b) =
        calc_scaling_factors(hap_a_kmers, hap_b_kmers);

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
        let read = seq::SeqRecord {
            id: "test".to_string(),
            seq: "ACGGGCATCGCGGC".to_string(),
            entry_string: ">test\nACGGGCATCGCGGC".to_string(),
        };

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
