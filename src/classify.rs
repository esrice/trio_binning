use kmer;
use seq;
use std::{result, error, fmt, cmp};
use std::fs::File;
use std::io::Write;

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

///
/// # Errors
/// io::Error -- if any input or output file can't be opened
/// seq::ExtensionError -- if the input file type cannot be determined
pub fn classify_unpaired(hap_a_kmers: &kmer::KmerSet,
                         hap_b_kmers: &kmer::KmerSet,
                         input_reads_filename: &str,
                         hap_a_out_prefix: &str,
                         hap_b_out_prefix: &str,
                         hap_u_out_prefix: &str,
                         k: usize) -> Result<()> {

    // set up input stream
    // this can return io::Error or seq::ExtensionError
    let input_reader = seq::SeqReader::from_path(input_reads_filename)?;

    // set up output streams
    // TODO figure out correct extensions based on input_reader type
    // TODO gzip output support
    let mut hap_a_out = File::create(hap_a_out_prefix)?; // io::Error
    let mut hap_b_out = File::create(hap_b_out_prefix)?; // io::Error
    let mut hap_u_out = File::create(hap_u_out_prefix)?; // io::Error

    // calculate read-count scaling factors
    let (scaling_factor_a, scaling_factor_b) =
        calc_scaling_factors(hap_a_kmers, hap_b_kmers);

    for result in input_reader {
        let record = result?;
        let (hap_a_count, hap_b_count) = count_kmers_in_read(
            hap_a_kmers, hap_b_kmers, &record, k)?;
        let hap_a_score = (hap_a_count as f32) * scaling_factor_a;
        let hap_b_score = (hap_b_count as f32) * scaling_factor_b;

        let mut haplotype = "?";
        if hap_a_score > hap_b_score {
            hap_a_out.write(record.entry_string.as_bytes())?;
            haplotype = "A";
        } else if hap_b_score > hap_a_score {
            hap_b_out.write(record.entry_string.as_bytes())?;
            haplotype = "B";
        } else {
            hap_u_out.write(record.entry_string.as_bytes())?;
            haplotype = "U";
        }
        println!("{}\t{}\t{}\t{}", record.id, haplotype,
                 hap_a_score, hap_b_score);
    }

    Ok(())
}

pub fn classify_paired(hap_a_kmers: &kmer::KmerSet, hap_b_kmers: &kmer::KmerSet,
                   input_reads_filename_a: &str, input_reads_filename_b: &str,
                   hap_a_out_prefix: &str, hap_b_output_prefix: &str,
                   hap_u_out_prefix: &str) -> Result<()> {

    // calculate read-count scaling factors
    let (scaling_factor_a, scaling_factor_b) =
        calc_scaling_factors(hap_a_kmers, hap_b_kmers);
    unimplemented!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count() {
        let read = seq::SeqRecord {
            id: "test".to_string(),
            seq: "AACAACGCGCGTCGGTATCT".to_string(),
            entry_string: ">test\nAACAACGCGCGTCGGTATCT".to_string(),
        };

        let hap_a_kmer_strings = ["AACAA", "AGATA", "TTTTT"];
        let hap_b_kmer_strings = ["GATTT", "AGTCG", "AACGC"];
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
